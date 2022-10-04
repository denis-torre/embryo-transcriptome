#################################################################
#################################################################
############### Embryo IsoSeq ################
#################################################################
#################################################################
##### Author: Denis Torre
##### Affiliation: Sebra and Guccione Laboratories,
##### Icahn School of Medicine at Mount Sinai

#############################################
########## 1. Load libraries
#############################################
##### 1. Python modules #####
from ruffus import *
from ruffus.combinatorics import *
import ruffus.cmdline as cmdline
import sys
import os
import json
import glob
import sqlite3
import pandas as pd
import numpy as np
# from rpy2.robjects import r, pandas2ri
# pandas2ri.activate()

##### 2. LSF #####
# 2.1 Import
sys.path.append('/hpc/users/torred23/pipelines/support')
import lsf

# 2.2 Default parameters
r_source = 'pipeline/scripts/embryo-isoseq.R'
py_source = 'pipeline/scripts/EmbryoIsoseq.py'
P = 'acc_apollo'
# q = 'express'
# q = 'premium'
q = 'sla'
W = '00:30'
GB = 5
n = 1
mkdir_val = True

# 2.3 Wrappers
# CMD
def run_job(cmd_str, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_job(cmd_str, outfile, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

# R
def run_r_job(func_name, func_input, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_r_job(func_name, func_input, outfile, r_source=r_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

# Py
def run_py_job(func_name, func_input, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_py_job(func_name, func_input, outfile, py_source=py_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

##### 3. Custom script imports #####
# 3.1 Python
#sys.path.append('pipeline/scripts')
#import EmbryoIsoseq as P

# 3.2 R
# r.source(r_source)

#############################################
########## 2. General Setup
#############################################
##### 1. Pipeline running #####
# Pipeline args
options = cmdline.get_argparse().parse_args()

##### 2. Genome indices #####
# Open JSON
with open('/hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/reference_genomes/reference-genomes.json') as openfile:
	genome_indices = json.load(openfile)


##### 3. Variables #####
# Primers
clontech_primer_file = 'arion/isoseq/s01-isoseq.dir/primer.fasta'

# Qiao
qiao_metadata = 'arion/datasets/qiao/qiao-sample_names.csv'
qiao_subreads = 'arion/datasets/qiao/rawdata/pacbio/*/*.subreads.bam'
qiao_flnc = 'arion/isoseq/s01-isoseq.dir/mouse/flnc/*.flnc.bam'

# Human
human_metadata = 'arion/datasets/human/human-sample_names.csv'
human_flnc = 'arion/datasets/human/human_embryo_pacb/*.flnc.bam'

# Illumina
illumina_fastq = 'arion/illumina/s01-fastq.dir/*/trimmed/*/*.fq.gz'

# Outlier samples
outlier_sample_file = '../embryo-illumina/pipeline/outlier_samples.json'
with open(outlier_sample_file) as openfile:
	outlier_samples = json.load(openfile)

#######################################################
#######################################################
########## S1. Process data
#######################################################
#######################################################

#############################################
########## 1. Index
#############################################

@transform(qiao_subreads,
		   suffix('.bam'),
		   '.bam.pbi')

def indexSubreads(infile, outfile):

	# Command
	cmd_str = ''' pbindex {infile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='isoseq3', W='00:15', GB=10, n=1, print_outfile=False, print_cmd=False)

#############################################
########## 2. CCS
#############################################

@follows(indexSubreads)

@subdivide(qiao_subreads,
		   regex(r'(.*)/datasets/(.*)/rawdata/pacbio/(.*)/.*.subreads.bam'),
		   r'\1/isoseq/s01-isoseq.dir/mouse/ccs/\3/chunk_*/\3_chunk_*.ccs.bam',
		   r'\1/isoseq/s01-isoseq.dir/mouse/ccs/\3/chunk_{i}/\3_chunk_{i}.ccs.bam')

def runCCS(infile, outfiles, outfileRoot):

	# Loop
	n = 50
	for i in range(1, n+1):

		# Get outfile
		outfile = outfileRoot.format(**locals())

		# Command
		cmd_str = 'ccs {infile} {outfile} --min-rq 0.9 --chunk {i}/{n}'.format(**locals())

		# Run
		run_job(cmd_str, outfile, conda_env='isoseq3', W='01:00', GB=4, n=6, print_outfile=False)

#############################################
########## 3. lima
#############################################
# primer.fasta file from https://github.com/PacificBiosciences/IsoSeq/blob/master/isoseq-clustering.md#step-2---primer-removal-and-demultiplexing on 2020/11/05

@transform(runCCS,
		   suffix('ccs.bam'),
		   add_inputs(clontech_primer_file),
		   'demux.Clontech_5p--NEB_Clontech_3p.bam')

def runLima(infiles, outfile):

	# Get basename
	basename = outfile.replace('.Clontech_5p--NEB_Clontech_3p', '')

	# Command
	cmd_str = 'lima --isoseq --different --min-passes 1 --split-named --dump-clips --dump-removed --peek-guess {infiles[0]} {infiles[1]} {basename}'.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='isoseq3', W='00:15', GB=10, n=1, print_cmd=False)

#############################################
########## 4. Refine
#############################################

@transform(runLima,
		   suffix('.bam'),
		   add_inputs(clontech_primer_file),
		   '.flnc.bam')

def refineReads(infiles, outfile):

	# Command
	cmd_str = 'isoseq3 refine --require-polya {infiles[0]} {infiles[1]} {outfile}'.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='isoseq3', W='00:15', GB=10, n=1, print_cmd=False)

#############################################
########## 5. Merge reads
#############################################

@collate(refineReads,
		 regex(r'(.*)/ccs/(.*)/chunk_.*/.*.flnc.bam'),
		 r'\1/flnc/\2.flnc.bam')

def mergeChunks(infiles, outfile):

	# Command
	infiles_str = ' '.join(infiles)
	cmd_str = 'pbmerge -o {outfile} {infiles_str}'.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='isoseq3', W='00:30', GB=15, n=1, print_outfile=False)

#############################################
########## 6. Merge reports
#############################################

# @follows(refineReads)

# @collate('arion/isoseq/s01-isoseq.dir/*/ccs/*/*/*.demux.Clontech_5p--NEB_Clontech_3p.flnc.report.csv',
# 		 regex(r'(.*)/ccs/(.*)/.*/.*.report.csv'),
# 		 r'\1/flnc/\2.flnc.report.csv')

# def mergeReports(infiles, outfile):

# 	# Read
# 	print('Doing {}...'.format(outfile))
# 	dataframe = pd.concat([pd.read_csv(x) for x in infiles])

# 	# Write
# 	dataframe.to_csv(outfile, index=False)

#############################################
########## 7. Convert to FASTA
#############################################

# @transform('arion/isoseq/s01-isoseq.dir/mouse/flnc/SRR10267008.flnc.bam',
@transform((mergeChunks, human_flnc),
		   suffix('.bam'),
		   '.fastq.gz')

def flncToFastq(infile, outfile):

	# Command
	cmd_str = ''' samtools view {infile} | awk '{{printf("@%s\\n%s\\n+\\n%s\\n", $1, $10, $11)}}' | gzip > {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['samtools/1.11'], W='01:00', GB=10, n=1, print_outfile=True, print_cmd=False)

#######################################################
#######################################################
########## S2. Align
#######################################################
#######################################################

#############################################
########## 1. Copy FASTQ
#############################################

def fastqJobs():

	# Paths
	paths = {
		'mouse': 'arion/isoseq/s01-isoseq.dir/mouse/flnc/{sample_id}.flnc.fastq.gz',
		'human': 'arion/datasets/human/human_embryo_pacb/{sample_id}.flnc.fastq.gz'
	}

	# Get sample names
	sample_dict = {x.split('/')[-2].replace('qiao', 'mouse'): pd.read_csv(x).query('Platform == "PACBIO_SMRT"').set_index('sample_id')['sample_name'].to_dict() for x in [qiao_metadata, human_metadata]}

	# Loop
	for organism, samples in sample_dict.items():
		for sample_id, sample_name in samples.items():
			infile = os.path.join(os.getcwd(), paths[organism].format(**locals()))
			outfile = 'arion/isoseq/s02-alignment.dir/fastq/{organism}/{sample_name}.flnc.fastq'.format(**locals())
			yield [infile, outfile]

@follows(flncToFastq)

@files(fastqJobs)

def copyFASTQ(infile, outfile):

	# Run
	run_r_job('copy_fastq', infile, outfile, run_locally=True)#conda_env='env', modules=[], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 2. minimap2
#############################################

@transform(copyFASTQ,
		   regex(r'(.*)/fastq/(.*)/(.*).flnc.fastq'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.primary_assembly.fa'),
		   r'\1/minimap2/\2/\3/\3-minimap2.sam')

def runMinimap2(infiles, outfile):

	# Command
	outname = outfile.rsplit('.', 1)[0]
	cmd_str = ''' minimap2 -ax splice -uf --secondary=no -C5 -t 30 --MD {infiles[1]} {infiles[0]} | samtools sort -O sam -o {outfile} && \
		samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | cut -f1,5 > {outname}.mapq '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['samtools/1.9', 'minimap2/2.17'], W='03:00', GB=6, n=6, print_outfile=False, stderr=outfile.replace('.sam', '.log'))

#############################################
########## 3. Filter
#############################################

@transform(runMinimap2,
		   regex(r'(.*).sam'),
		   r'\1.filtered.sam')

def filterSam(infile, outfile):

	# Files
	outname = outfile.replace('.sam', '')

	# Command
	cmd_str = ''' sambamba view --with-header --nthreads 30 --sam-input --format sam --filter "not unmapped and mapping_quality >= 50" {infile} | samtools sort -O sam -o {outfile} && \
		samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | cut -f5 -d "	" > {outname}.mapq'''.format(**locals())
	# and not duplicate

	# Run
	run_job(cmd_str, outfile, W='00:15', n=5, GB=1, modules=['samtools/1.9', 'sambamba/0.5.6'], print_cmd=False)

#############################################
########## 4. FastQC
#############################################

@transform('arion/isoseq/s02-alignment.dir/fastq/*/*.flnc.fastq',
		   regex(r'(.*)/fastq/(.*)/(.*).fastq'),
		   r'\1/fastqc/\2/\3.fastqc.html')

def runFastQC(infile, outfile):

	# Command
	cmd_str = '''fastqc --outdir=$(dirname {outfile}) {infile}'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['fastqc/0.11.8'], W='03:00', GB=12, n=1, print_outfile=False)

#############################################
########## 5. MultiQC
#############################################

@transform('arion/isoseq/s02-alignment.dir/fastqc/*',
		   regex(r'(.*)/.*.dir/fastqc/(.*)'),
		   r'\1/multiqc/fastqc/\2/multiqc_report.html')

def runMultiQC(infile, outfile):

	# Command
	cmd_str = 'multiqc --outdir $(dirname {outfile}) {infile}'.format(**locals())

	# Run
	if not os.path.exists(outfile):
		run_job(cmd_str, outfile, conda_env='env', W="00:20", GB=5, n=1, print_outfile=False, run_locally=False)

#######################################################
#######################################################
########## S3. Illumina alignment
#######################################################
#######################################################

#############################################
########## 1. STAR index
#############################################

@follows(filterSam)

@transform('arion/datasets/reference_genomes/*/*.102.gtf',
		   regex(r'(.*)/(.*)/(.*).102.gtf'),
		   add_inputs(r'\1/\2/\3.dna_sm.primary_assembly.fa'),
		   r'arion/isoseq/s03-illumina_alignment.dir/\2/STAR/index/')

def buildStarIndex(infiles, outfile):

	# Split
	reference_gtf, reference_fasta = infiles

	# Command
	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {reference_fasta} --sjdbGTFfile {reference_gtf} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='01:00', GB=5, n=15, ow=True, print_cmd=False, jobname='_'.join(outfile.split('/')[-4:-1]))

#############################################
########## 2. STAR pass 1
#############################################

@follows(buildStarIndex)

@collate(illumina_fastq,
		regex(r'.*/s01-fastq.dir/(.*)/trimmed/(.*)/.*.fq.gz'),
		add_inputs(r'arion/isoseq/s03-illumina_alignment.dir/\1/STAR/index/'),
		r'arion/isoseq/s03-illumina_alignment.dir/\1/STAR/pass1/\2/\2-SJ.out.tab')

def getSJsPass1(infiles, outfile):

	# Prefix
	prefix = outfile[:-len('SJ.out.tab')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {infiles[0][1]} \
		--readFilesIn {infiles[0][0]} {infiles[1][0]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 50 \
		--outSAMtype None'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="02:00", GB=5, n=10, modules=['star/2.7.5b'], q='express')

# ls arion/isoseq/s03-illumina_alignment.dir/*/STAR/pass1/*/*-Log.progress.out | lr

#############################################
########## 3. STAR pass 2
#############################################

@follows(getSJsPass1)

@collate(illumina_fastq,
		regex(r'.*/s01-fastq.dir/(.*)/trimmed/(.*)/.*.fq.gz'),
		add_inputs(r'arion/isoseq/s03-illumina_alignment.dir/\1/STAR/index/', r'arion/isoseq/s03-illumina_alignment.dir/\1/STAR/pass1/*/*-SJ.out.tab'),
		r'arion/isoseq/s03-illumina_alignment.dir/\1/STAR/pass2/\2/\2-Aligned.sortedByCoord.out.bam')

def getSJsPass2(infiles, outfile):

	# Prefix
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
	sj_files_str = ' '.join(infiles[0][2:])

	# Command
	cmd_str = ''' STAR \
		--genomeDir {infiles[0][1]} \
		--readFilesIn {infiles[0][0]} {infiles[1][0]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 50 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 5000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="10:00", GB=10, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_cmd=True, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.log'))

#############################################
########## 4. Merge SJ files
#############################################

# @follows(getSJsPass1, getSJsPass2)

@collate(('arion/isoseq/s03-illumina_alignment.dir/*/STAR/pass1/*/*-SJ.out.tab', 'arion/isoseq/s03-illumina_alignment.dir/*/STAR/pass2/*/*-SJ.out.tab'),
		 regex(r'(.*.dir)/(.*)/STAR/.*.tab'),
		 r'\1/\2/STAR/\2-merged.SJ.out.tab')

def mergeSJs(infiles, outfile):

	# Run
	run_r_job('merge_sjs', infiles, outfile, conda_env='env', W='00:15', GB=15, n=1, stdout=outfile.replace('.tab', '.log'), stderr=outfile.replace('.tab', '.err'))#, q='premium')

#######################################################
#######################################################
########## S4. TranscriptClean
#######################################################
#######################################################

#############################################
########## 1. Trim genome FASTA
#############################################

@transform('arion/datasets/reference_genomes/*/*.dna_sm.primary_assembly.fa',
		   suffix('.fa'),
		   '_renamed.fa')

def renameFastaChromosomes(infile, outfile):

	# Command
	cmd_str = ''' cut -d ' ' -f1 {infile} > {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, W="00:30", GB=15, n=1)

#############################################
########## 2. Clean
#############################################

# @follows(mergeSJs, renameFastaChromosomes)

@transform(filterSam,
		   regex(r'(.*)/s02-alignment.dir/minimap2/(.*)/(.*)/.*.sam'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.dna_sm.primary_assembly_renamed.fa', r'arion/isoseq/s03-illumina_alignment.dir/\2/STAR/\2-merged.SJ.out.tab'),
		   r'\1/s04-cleaned_transcripts.dir/\2/\3/\3_clean.sam')

def runTranscriptClean(infiles, outfile):

	# Prefix
	tempdir = os.path.join(os.path.dirname(outfile), 'tmp')
	prefix = outfile[:-len('_clean.sam')]

	# Command
	cmd_str = ''' python /sc/arion/work/torred23/libraries/TranscriptClean/TranscriptClean-2.0.2/TranscriptClean.py \
		--sam {infiles[0]} \
		--genome {infiles[1]} \
		--spliceJns {infiles[2]} \
		--tmpDir {tempdir} \
		--threads 10 \
		--deleteTmp \
		--canonOnly \
		--primaryOnly \
		--outprefix {prefix}
	'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, print_cmd=False, W="03:00", GB=5, n=10, conda_env='talon', stdout=outfile.replace('_clean.sam', '_job.log'), stderr=outfile.replace('_clean.sam', '_job.err'))

#############################################
########## 3. Flag reads
#############################################

@transform(runTranscriptClean,
		   regex(r'(.*.dir)/(.*)/(.*)/(.*).sam'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.dna_sm.primary_assembly_renamed.fa'),
		   r'\1/\2/\3/flagged/\4_labeled.sam')

def flagReads(infiles, outfile):

	# Prefix
	tempdir = os.path.join(os.path.dirname(outfile), 'tmp')
	prefix = outfile[:-len('_labeled.sam')]

	# Command
	cmd_str = ''' talon_label_reads \
		--f {infiles[0]} \
		--g {infiles[1]} \
		--t 1 \
		--ar 20 \
		--tmpDir {tempdir} \
		--deleteTmp \
		--o {prefix}
	'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, W="02:00", GB=15, n=1, conda_env='talon', stdout=outfile.replace('_labeled.sam', '_job.log'), stderr=outfile.replace('_labeled.sam', '_job.err'))

#######################################################
#######################################################
########## S5. TALON
#######################################################
#######################################################

#############################################
########## 1. Initialize
#############################################

@transform('arion/datasets/reference_genomes/*/*GRC*.102.gtf',
		   regex(r'(.*)/(.*)/(.*)(.102).gtf'),
		   add_inputs(r'\1/\2/\3.dna_sm.primary_assembly.fa'),
		   r'arion/isoseq/s05-talon.dir/\2/\3\4_talon.db')

def initializeTalonDatabase(infiles, outfile):

	# Get parameters
	annotation_name = os.path.basename(outfile)[:-len('.db')]
	genome_build = annotation_name.split('.')[1].replace('GRCm38', 'mm10').replace('GRCh38', 'hg38')
	outname = outfile[:-len('.db')]

	# Command
	cmd_str = ''' talon_initialize_database \
		--f {infiles[0]} \
		--a {annotation_name} \
		--g {genome_build} \
		--l 200 \
		--idprefix TALON \
		--5p 1000 \
		--3p 1000 \
		--o {outname} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='01:00', n=1, GB=15, conda_env='talon', print_cmd=False, stdout=outfile.replace('.db', '.log'), stderr=outfile.replace('.db', '.err'), wait=False)

#############################################
########## 2. Config
#############################################

@collate(flagReads,
		 regex(r'(.*)/s04-cleaned_transcripts.dir/(.*?)/.*.sam'),
		 r'\1/s05-talon.dir/\2/\2-talon_config.csv')

def createTalonConfig(infiles, outfile):

	# Create dataframe
	config_dataframe = pd.DataFrame([{
		'dataset_name': x.split('/')[-3],
		'sample_description': x.split('/')[-3],
		'platform': 'PacBio-Sequel2' if 'human' in outfile else 'PacBio-Sequel',
		'sam_file': os.path.join(os.getcwd(), x),
	} for x in infiles])

	# outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Write
	config_dataframe.to_csv(outfile, header=False, index=False)

#############################################
########## 3. Run
#############################################

# @follows(createTalonConfig)

@transform(initializeTalonDatabase,
		   regex(r'(.*)/(.*).db'),
		   add_inputs(r'\1/*-talon_config.csv'),
		   r'\1/\2_annotated_QC.log')

def runTalon(infiles, outfile):

	# Get parameters
	annotation_name = os.path.basename(infiles[0])[:-len('.db')]
	genome_build = annotation_name.split('.')[1].replace('GRCm38', 'mm10').replace('GRCh38', 'hg38')

	# Fix paths
	outdir = os.path.dirname(outfile)+'/'
	infiles = [os.path.basename(x) for x in infiles]
	outname = os.path.basename(outfile[:-len('_QC.log')])

	# Command
	cmd_str = ''' cd {outdir} && talon \
		--db {infiles[0]} \
		--f {infiles[1]} \
		--build {genome_build} \
		--threads 10 \
		--cov 0.99 \
		--identity 0.95 \
		--o {outname} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='06:00', n=10, GB=10, conda_env='talon', print_cmd=True, stdout=outfile.replace('_QC.log', '.log'), stderr=outfile.replace('_QC.log', '.err'))

#############################################
########## 4. Filter transcripts
#############################################

# @follows(runTalon)

@transform(initializeTalonDatabase,
		   regex(r'(.*).db'),
		   r'\1_transcripts.tsv')

def filterTranscripts(infile, outfile):
	
	# Open database
	with sqlite3.connect(infile) as conn:

		# Query
		query = """ SELECT DISTINCT gene_ID, transcript_ID 
						FROM transcripts
						WHERE transcript_ID NOT IN (
							SELECT DISTINCT transcript_ID
								FROM observed
								GROUP BY transcript_ID
								HAVING MIN(fraction_As) > 0.6
						)
		""".format(**locals())

		# Get dataframe
		id_dataframe = pd.read_sql_query(query, conn)

	# Write
	id_dataframe.to_csv(outfile, header=False, index=False)

#############################################
########## 5. Get GTF
#############################################

@follows(filterTranscripts)

@transform(initializeTalonDatabase,
		   regex(r'(.*)/(.*).db'),
		   add_inputs(r'\1/\2_transcripts.tsv'),
		   r'\1/gtf/\2.gtf')

def getTalonGTF(infiles, outfile):
	
	# Get parameters
	annotation_name = os.path.basename(infiles[0])[:-len('.db')]
	genome_build = annotation_name.split('.')[1].replace('GRCm38', 'mm10').replace('GRCh38', 'hg38')
	prefix = outfile[:-len('_talon.gtf')]

	# Command
	cmd_str = ''' talon_create_GTF \
		--db {infiles[0]} \
		--annot {annotation_name} \
		--build {genome_build} \
		--whitelist {infiles[1]} \
		--o {prefix} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='00:10', n=1, GB=10, conda_env='talon', print_cmd=False, stdout=outfile.replace('.gtf', '_gtf.log'), stderr=outfile.replace('.gtf', '_gtf.err'))

#############################################
########## 6. Get abundance
#############################################

# @follows(filterTranscripts)

@transform('arion/isoseq/s05-talon.dir/*/*.102_talon.db',
		   regex(r'(.*)/(.*).db'),
		   add_inputs(r'\1/\2_transcripts.tsv'),
		   r'\1/\2_abundance_filtered.tsv')

def getTalonAbundance(infiles, outfile):
	
	# Get parameters
	annotation_name = os.path.basename(infiles[0])[:-len('.db')]
	genome_build = annotation_name.split('.')[1].replace('GRCm38', 'mm10').replace('GRCh38', 'hg38')
	prefix = outfile[:-len('_talon_abundance_filtered.tsv')]

	# Command
	cmd_str = ''' talon_abundance \
		--db {infiles[0]} \
		--annot {annotation_name} \
		--build {genome_build} \
		--whitelist {infiles[1]} \
		--o {prefix} '''.format(**locals())
		# --observed \

	# Run
	run_job(cmd_str, outfile, W='01:00', n=1, GB=25, conda_env='talon', print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 7. SQANTI
#############################################

@transform(getTalonGTF,
		   regex(r'(.*)/(.*)/gtf/(.*).gtf'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.gtf', r'arion/datasets/reference_genomes/\2/*.primary_assembly.fa', r'/sc/arion/projects/GuccioneLab/genome-indices/*/polya/\2_polyA_motifs.txt', r'/sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.2/data/ref_TSS_annotation/\2.refTSS_v3.1.????.nochr.bed', r'arion/isoseq/s03-illumina_alignment.dir/\2/STAR/pass2/*/*-SJ.out.tab'),
		   r'\1/\2/gtf/sqanti/\3_SQANTI3_report.pdf')

def runSqanti(infiles, outfile):
	
	# Add paths
	transcript_gtf, ensembl_gtf, ensembl_fasta, polya_file, tss_file = infiles[:5]
	sj_files = ','.join([x for x in infiles[5:] if not any([y in x for y in outlier_samples[outfile.split('/')[-4]]])])

	# Get output
	outdir = os.path.dirname(outfile)
	output = os.path.basename(outfile)[:-len('_SQANTI3_report.pdf')]

	# PYTHONPATH
	PYTHONPATH = '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.2/cDNA_Cupcake/:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.2/cDNA_Cupcake/sequence'

	# Command
		# --cage_peak {tss_file} \
	cmd_str = ''' export PYTHONPATH={PYTHONPATH} && /sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.2/sqanti3_qc.py {transcript_gtf} {ensembl_gtf} {ensembl_fasta} \
		--polyA_motif_list {polya_file} \
		--coverage {sj_files} \
		--force_id_ignore \
        --skipORF \
        --report pdf \
		--dir {outdir} \
		--output {output}
	'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, conda_env='SQANTI3_v4.2', W='06:00', GB=50, n=1, print_cmd=False, stderr=outfile.replace('.pdf', '.err'), stdout=outfile.replace('.pdf', '.log'))

#############################################
########## 7. Get SJ support
#############################################

@follows(runSqanti)

@transform('arion/isoseq/s05-talon.dir/*/gtf/sqanti/*.102_talon_junctions.txt',
		   suffix('s.txt'),
		   '_counts.tsv')

def getSjCounts(infile, outfile):

	# Run
	run_r_job('get_sj_counts', infile, outfile, conda_env='env', W='06:00', GB=30, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 8. Filter GTF
#############################################

@follows(getTalonAbundance, getSjCounts)

@transform(getTalonGTF,
		   regex(r'(.*)/(.*)/(.*).gtf'),
		   add_inputs(r'\1/\2/sqanti/\3_junction_counts.tsv', r'\1/\3_abundance_filtered.tsv'),
		   r'\1/\2/\3-SJ_filtered.gtf')

def filterGTF(infiles, outfile):

	# Run
	run_r_job('filter_gtf', infiles, outfile, print_outfile=False, conda_env='env', W='00:30', GB=10, n=1, stdout=outfile.replace('.gtf', '.log'), stderr=outfile.replace('.gtf', '.err'))#, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 9. Get FASTA
#############################################

@transform(filterGTF,
		   regex(r'(.*)/(.*)/(.*)/(.*).gtf'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.dna_sm.primary_assembly.fa'),
		   r'\1/\2/\3/\4.fasta')

def getFASTA(infiles, outfile):

	# Command
	cmd_str = ''' gffread {infiles[0]} -g {infiles[1]} -w {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['gff/2021-02'], W='00:30', GB=10, n=1, stdout=outfile.replace('.fasta', '_fasta.log'), stderr=outfile.replace('.fasta', '_fasta.err'))

#############################################
########## 10. Classify transcripts
#############################################

@follows(getTalonAbundance, runSqanti, filterGTF)

@transform('arion/isoseq/s05-talon.dir/*/gtf/*.gtf',
		   regex(r'(.*)/gtf/(.*talon)(.*).gtf'),
		   add_inputs(r'\1/\2_abundance_filtered.tsv', r'\1/gtf/sqanti/\2_classification.txt'),
		   r'\1/gtf/\2\3-transcript_classification.tsv')

def createTranscriptClassification(infiles, outfile):

	# Run
	run_r_job('create_transcript_classification', infiles, outfile, conda_env='env', W='00:15', GB=10, n=1, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 10. Split transcript GTF
#############################################

@follows(getTalonGTF, filterGTF)

@subdivide('arion/isoseq/s06-cpat.dir/human/gtf/*-SJ_filtered.cds.gtf',
		   regex(r'(.*)/gtf/(.*talon)(.*).gtf'),
		   add_inputs(r'arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered-transcript_classification.tsv'),
		   r'\1/gtf/split_classification\3/\2\3_*.gtf',
		   r'\1/gtf/split_classification\3/\2\3_{transcript_class}.gtf')

def splitTranscriptGTF(infiles, outfiles, outfileRoot):

	# Run
	run_r_job('split_transcript_gtf', infiles, outfileRoot, conda_env='env', W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 11. Sort and index
#############################################

@follows(splitTranscriptGTF)

@transform('arion/isoseq/s06-cpat.dir/human/gtf/split_classification-SJ_filtered.cds/Homo_sapiens.GRCh38.102_talon-SJ_filtered.cds_*.gtf',
		   suffix('.gtf'),
		   '.sorted.gtf')

def indexGTF(infile, outfile):

	# Command
	cmd_str = ''' igvtools sort -m 5000000 {infile} {outfile} && igvtools index {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env=False, modules=['igvtools/2.3.32'], W='00:15', GB=20, n=1, print_outfile=False, print_cmd=False)

#############################################
########## 11. Get junctions
#############################################

@transform('arion/isoseq/s05-talon.dir/human/gtf/*_talon-SJ_filtered.gtf',
		   suffix('.gtf'),
		   '_junctions.tsv')

def getJunctions(infile, outfile):

	# Run
	run_r_job('get_junctions', infile, outfile, conda_env='env', W='03:00', GB=10, n=5, run_locally=True, print_outfile=False, print_cmd=False)

#############################################
########## 11. Annotate splice sites
#############################################

@transform('arion/isoseq/s05-talon.dir/human/gtf/*_talon-SJ_filtered.gtf',
		   suffix('.gtf'),
		   '_splice_sites.tsv')

def annotateSpliceSites(infile, outfile):

	# Run
	run_r_job('annotate_splice_sites', infile, outfile, conda_env='env', W='03:00', GB=30, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#######################################################
#######################################################
########## S6. CPAT
#######################################################
#######################################################

#############################################
########## 1. Run
#############################################

@follows(getFASTA)

@transform('arion/isoseq/s05-talon.dir/*/gtf/*_talon-SJ_filtered.fasta',
		   regex(r'(.*)/s05-talon.dir/(.)(.*)/gtf/(.*).fasta'),
		   add_inputs(r'/sc/arion/work/torred23/libraries/CPAT3/CPAT-3.0.0/dat/*\3_Hexamer.tsv', r'/sc/arion/work/torred23/libraries/CPAT3/CPAT-3.0.0/dat/*\3_logitModel.RData'),
		   r'\1/s06-cpat.dir/\2\3/\4-cpat.ORF_prob.best.tsv')

def runCPAT(infiles, outfile):
	
	# Split
	transcript_fasta, cpat_hexamer, cpat_model = infiles
	basename = outfile.replace('.ORF_prob.best.tsv', '')

	# Command
	cmd_str = ''' ~/.conda/envs/env/bin/cpat.py -g {transcript_fasta} \
		-x {cpat_hexamer} \
		-d {cpat_model} \
		--top-orf=5 \
		--log-file={basename}.log \
		-o {basename}
	'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, conda_env='env', W='10:00', GB=25, n=1, print_outfile=True, wait=False)

#############################################
########## 2. Get results
#############################################

@transform(runCPAT,
		   suffix('.tsv'),
		   '_results.tsv')

def formatCPAT(infile, outfile):

	# Rename
	rename_dict = {'seq_ID': 'Sequence Name', 'mRNA': 'RNA size', 'ORF': 'ORF size', 'Fickett': 'Ficket Score', 'Hexamer': 'Hexamer Score', 'Coding_prob': 'Coding Probability'}

	# Coding cutoff (https://cpat.readthedocs.io/en/latest/ on 2021/03/17)
	coding_cutoff = 0.364 if 'human' in outfile else 0.44

	# Read
	cpat_dataframe = pd.read_table(infile).rename(columns=rename_dict)[rename_dict.values()]
	cpat_dataframe['Coding Label'] = ['yes' if x >= coding_cutoff else 'no' for x in cpat_dataframe['Coding Probability']]
	cpat_dataframe.index.name = 'Data ID'

	# Write
	cpat_dataframe.to_csv(outfile, sep='\t', index=True)

#############################################
########## 3. Split GTF
#############################################

@subdivide(filterGTF,
		   regex(r'(.*)/s05-talon.dir/(.*)/gtf/(.*).gtf'),
		   r'\1/s06-cpat.dir/\2/gtf/split/\3_??.gtf',
		   r'\1/s06-cpat.dir/\2/gtf/split/\3_{chunk_nr}.gtf')

def splitGTF(infile, outfiles, outfileRoot):

	# Run
	if not len(outfiles):
		run_r_job('split_gtf', infile, outfileRoot, conda_env='env', W='01:00', GB=15, n=1, run_locally=False, wait=False)

#############################################
########## 4. Add CDS
#############################################

@follows(formatCPAT)

@transform(splitGTF,
		   regex(r'(.*)/(gtf/split/.*).gtf'),
		   add_inputs(r'\1/*-cpat.ORF_prob.best.tsv'),
		   r'\1/\2.cds.gtf')

def addCDS(infiles, outfile):

	# Coding cutoff (https://cpat.readthedocs.io/en/latest/ on 2021/03/17)
	coding_cutoff = 0.364 if 'human' in outfile else 0.44

	# Run
	run_r_job('add_cds', infiles, outfile, additional_params=coding_cutoff, conda_env='env', W='02:00', GB=25, n=1, run_locally=False, stdout=outfile.replace('.gtf', '.log'), stderr=outfile.replace('.gtf', '.err'))

#############################################
########## 5. Merge
#############################################

@collate(addCDS,
		 regex(r'(.*)/split/(.*)_.*.cds.gtf'),
		 r'\1/\2.cds.gtf')

def mergeGTF(infiles, outfile):

	# Run
	run_r_job('merge_gtf', infiles, outfile, conda_env='env', W='00:45', GB=15, n=1, run_locally=False)

#######################################################
#######################################################
########## S10. Pfam
#######################################################
#######################################################

#############################################
########## 1. Translate
#############################################

@transform(runCPAT,
		   regex(r'(.*)/(s06-cpat.dir)/(.*)/(.*).ORF_prob.best.tsv'),
		   add_inputs(r'\1/\2/\3/\4.ORF_seqs.fa'),
		   r'\1/s07-pfam.dir/\3/fasta/\3-translated.fasta')

def translateORFs(infiles, outfile):

	# Run
	run_r_job('translate_orfs', infiles, outfile, conda_env='env', W="01:00", GB=15, n=1, stdout=outfile.replace('.fasta', '.log'), wait=False)

#############################################
########## 2. Split
#############################################

@subdivide(translateORFs,
		   regex(r'(.*).fasta'),
		   r'\1.fasta.*',
		   r'\1.fasta.100')

def splitORFs(infile, outfiles, outfileRoot):

	# Get number
	N = outfileRoot.split('.')[-1]

	# Command
	cmd_str = ''' gt splitfasta -numfiles {N} {infile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfileRoot, modules=['genometools/1.5.9'], W='00:15', GB=10, n=1, jobname='splitORF_'+outfileRoot.split('/')[-3], wait=False)

#############################################
########## 3. Run
#############################################

@transform(splitORFs,
		   regex(r'(.*)/fasta/(.*).fasta\.(.*)'),
		   add_inputs('/sc/arion/work/torred23/libraries/PfamScan'),
		   r'\1/split/\2_\3_pfam.txt')

def runPfamScan(infiles, outfile):

	# Data directory
	input_fasta, pfam_dir = infiles

	# Command
	cmd_str = ''' {pfam_dir}/pfam_scan.pl \
		-dir {pfam_dir}/data \
		-cpu 50 \
		-fasta {input_fasta} > {outfile}
	'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, conda_env='env', modules=['hmmer/3.3'], W='06:00', GB=2, n=10, print_cmd=False, run_locally=False, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'))

#############################################
########## 3. Format
#############################################

@collate(runPfamScan,
		 regex(r'(.*)/split/(.*)_.*_pfam.txt'),
		 r'\1/\2_pfam.tsv')

def mergePfamResults(infiles, outfile):

	# Initialize
	results = []

	# Loop
	for infile in infiles:

		# Read
		pfam_dataframe = pd.read_csv(infile, comment='#', delim_whitespace=True, header=None)

		# Get column names
		colnames = [x.replace(' ', '_').replace('-', '_') for x in pd.read_csv(infile).iloc[26][0][:-1].replace('# ', '').replace('<', '').split('> ')]

		# Add column names
		pfam_dataframe.columns = colnames

		# Fix sequence ID
		pfam_dataframe['seq_id'] = [x.split('_ORF')[0] for x in pfam_dataframe['seq_id']]

		# Append
		results.append(pfam_dataframe)

	# Concatenate
	result_dataframe = pd.concat(results).query('E_value < 0.1').sort_values('seq_id')

	# Write
	result_dataframe.to_csv(outfile, index=False, sep='\t')
	
#######################################################
#######################################################
########## S8. RepeatMasker
#######################################################
#######################################################

#############################################
########## 1. Split
#############################################

@follows(getFASTA)

@subdivide('arion/isoseq/s05-talon.dir/*/gtf/*_talon-SJ_filtered.fasta',
		   regex(r'(.*)/s05-talon.dir/(.*)/gtf/(.*).fasta'),
		   r'\1/s08-repeatmasker.dir/\2/fasta/*.fasta.*',
		   r'\1/s08-repeatmasker.dir/\2/fasta/\3.fasta.100')

def splitFASTA(infile, outfiles, outfileRoot):

	# Get number
	N = outfileRoot.split('.')[-1]

	# Get temp filename
	outdir = os.path.dirname(outfileRoot)
	tempfile = os.path.join(outdir, os.path.basename(infile))

	# Command
	cmd_str = ''' cp {infile} {outdir} && gt splitfasta -numfiles {N} {tempfile} && rm {tempfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfileRoot, modules=['genometools/1.5.9'], W='00:15', GB=10, n=1, jobname='splitFASTA_'+outfileRoot.split('/')[-3], stdout=os.path.join(os.path.dirname(outfileRoot), 'job.log'))

#############################################
########## 2. Run
#############################################

@transform(splitFASTA,
		   regex(r'(.*)/fasta/(.*)'),
		   r'\1/split/\2.out')

def runRepeatMasker(infile, outfile):

	# Paths
	infile = os.path.join(os.getcwd(), infile)
	outdir = os.path.dirname(outfile)

	# Species
	species = 'Homo sapiens' if 'human' in outfile else 'Mus musculus'

	# Command
	cmd_str = ''' cd {outdir} && RepeatMasker -species "{species}" -pa 64 -dir . {infile}'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['repeatmasker/4.1.1'], W='06:00', GB=3, n=10, print_outfile=False, stdout=outfile.replace('.out', '.log'), stderr=outfile.replace('.out', '.err'))

#############################################
########## 3. Filter
#############################################

@collate(runRepeatMasker,
		 regex(r'(.*)/split/(.*talon).*.out'),
		 r'\1/\2_repeatmasker.tsv')

def mergeRepeatMasker(infiles, outfile):

	# Run
	run_r_job('merge_repeatmasker', infiles, outfile, conda_env='env', run_locally=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#######################################################
#######################################################
########## S9. Evolutionary conservation
#######################################################
#######################################################

#############################################
########## 1. Get BED
#############################################

# @follows(splitGTF)

@subdivide('arion/isoseq/s06-cpat.dir/human/gtf/split/*talon-SJ_filtered_??.gtf',
		   regex(r'(.*)/s06-cpat.dir/(.*)/gtf/split/(.*).gtf'),
		   r'\1/s09-evolutionary_conservation.dir/\2/bed/*/*.bed',
		   r'\1/s09-evolutionary_conservation.dir/\2/bed/{feature_type}/\3-{feature_type}.bed')

def gtfToBed(infile, outfiles, outfileRoot):

	# Loop
	for feature_type in ['exon', 'transcript']:

		# Get outfile
		outfile = outfileRoot.format(**locals())

		# Run
		run_r_job('gtf_to_bed', infile, outfile, additional_params=feature_type, conda_env='env', W='00:05', GB=10, n=1, stdout=outfile.replace('.bed', '.log'), stderr=outfile.replace('.bed', '.err'))

#############################################
########## 2. Get chromosome sizes
#############################################

@transform(renameFastaChromosomes,
		   suffix('.fa'),
		   '.chromsizes')

def getChromosomeSizes(infile, outfile):

	# Command
	cmd_str = ''' faidx {infile} -i chromsizes > {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['python/3.8.2'], W='00:05', GB=10, n=11, print_outfile=False)

#############################################
########## 3. Get excluded regions
#############################################

@transform('arion/datasets/reference_genomes/*/*.102.gtf',
		   regex(r'(.*).(...).gtf'),
		   add_inputs(r'\1*dna_sm.primary_assembly_gaps.bed'),
		   r'\1.\2_to_exclude_intergenic.bed')

def getExcludedRegions(infiles, outfile):

	# Run
	run_r_job('get_excluded_regions', infiles, outfile, conda_env='env', W='00:10', GB=10, n=1)#, stdout=outfile.replace('.bed', '.log'), stderr=outfile.replace('.bed', '.err'))

#############################################
########## 4. Shuffle transcripts
#############################################

@follows(gtfToBed, getChromosomeSizes, getKnownTranscriptBed)

@subdivide('arion/isoseq/s09-evolutionary_conservation.dir/*/bed/transcript/*.bed',
		   regex(r'(.*)/(.*)/bed/transcript/(.*).bed'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*renamed.chromsizes', r'arion/datasets/reference_genomes/\2/*to_exclude_intergenic.bed'),
		   r'\1/\2/bed_shuffled/*/transcript/\3-shuffled*.bed',
		   r'\1/\2/bed_shuffled/{shuffle_type}/transcript/\3-shuffled{shift_nr}-{shuffle_type}.bed')

def shuffleTranscripts(infiles, outfiles, outfileRoot):

	# Shuffle parameters
	shuffle_types = {
		'intergenic': '-excl {infiles[2]}'.format(**locals()),
		'intergenic_chrom': '-chrom -excl {infiles[2]}'.format(**locals())
	}

	# Loop
	for shuffle_type, shuffle_parameters in shuffle_types.items():

		# Loop
		for i in range(1):

			# Get shift number
			shift_nr = i+1

			# Get outfile
			outfile = outfileRoot.format(**locals())

			# Command
			cmd_str = ''' bedtools shuffle {shuffle_parameters} -i {infiles[0]} -g {infiles[1]} > {outfile} \n'''.format(**locals()) # -chrom
			
			# Run
			run_job(cmd_str, outfile, modules=['bedtools/2.29.2'], W='00:05', GB=5, n=1, print_cmd=False)

#############################################
########## 5. Get shuffled exons
#############################################

@transform(shuffleTranscripts,
		   regex(r'(.*)/bed_shuffled/(.*)/transcript/(.*)-transcript-(shuffled.*).bed'),
		   add_inputs(r'\1/bed/exon/\3-exon.bed'),
		   r'\1/bed_shuffled/\2/exon/\3-exon-\4.bed')

def getShuffledExons(infiles, outfile):

	# Run
	run_r_job('get_shuffled_exons', infiles, outfile, conda_env='env', W='00:05', GB=10, n=1)#, stdout=outfile.replace('.bed', '.log'), stderr=outfile.replace('.bed', '.err'))

#############################################
########## 6. Get scores
#############################################

def scoreJobs():
	for organism in ['human']: #mouse
		transcript_bed = glob.glob('arion/isoseq/s09-evolutionary_conservation.dir/{organism}/bed/exon/*.bed'.format(**locals()))
		shuffled_bed = glob.glob('arion/isoseq/s09-evolutionary_conservation.dir/{organism}/bed_shuffled/*/exon/*.bed'.format(**locals()))
		for score_file in glob.glob('arion/datasets/evolutionary_conservation/{organism}/*.bw'.format(**locals())):
			score_name = os.path.basename(score_file)[:-len('.bw')]
			for bed_file in transcript_bed+shuffled_bed:#: #
				bed_name = os.path.basename(bed_file)[:-len('.bed')]
				outfile = 'arion/isoseq/s09-evolutionary_conservation.dir/{organism}/scores/split/{score_name}/{bed_name}-{score_name}.tsv'.format(**locals())
				yield [[score_file, bed_file], outfile]

@files(scoreJobs)

def getConservationScores(infiles, outfile):

	# Command
	cmd_str = ''' bigWigAverageOverBed {infiles[0]} {infiles[1]} {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['ucsc-utils/2020-03-17'], W='00:05', GB=1, n=15, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 7. Merge
#############################################

@collate(getConservationScores,
		 regex(r'(.*)/scores/split/(.*)/.*talon-SJ_filtered_..(.*)-.*'),
		 r'\1/scores/average/\2/\2\3-average_scores.tsv')

def mergeConservationScores(infiles, outfile):

	# Run
	run_r_job('merge_conservation_scores', infiles, outfile, conda_env='env', W='00:10', GB=10, n=1, ow=True)

#############################################
########## 8. Get SS scores
#############################################

@transform('arion/datasets/evolutionary_conservation/human/*.bw',
		   regex(r'.*/(.*)/(.*).bw'),
		   add_inputs('arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered_splice_site_ranges_10bp.bed'),
		   r'arion/isoseq/s09-evolutionary_conservation.dir/\1/splice_site_conservation/\2_splice_sites_10bp.tsv')

def getSpliceSiteConservation(infiles, outfile):

	# Command
	cmd_str = ''' bigWigAverageOverBed {infiles[0]} {infiles[1]} {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['ucsc-utils/2020-03-17'], W='00:05', GB=1, n=15, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#######################################################
#######################################################
########## S10. liftOver
#######################################################
#######################################################

#############################################
########## 1. Convert
#############################################

@follows(addCDS)

@transform('arion/isoseq/s06-cpat.dir/*/gtf/split/*.cds.gtf',
		   regex(r'(.*)/s06-cpat.dir/(.*)/gtf/split/(.*).gtf'),
		   r'\1/s10-liftover.dir/\2/gp/\3.gp')

def convertToGenePred(infile, outfile):

	# Command
	tempfile = outfile.replace('.gp', '.tmp.gp')
	cmd_str = ''' ldHgGene -gtf -nobin -out={tempfile} ignored "" {infile} && awk 'OFS="\t" {{$2="chr"$2; print}}' {tempfile} > {outfile} && rm {tempfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['ucsc-utils/2020-03-17'], W='00:05', GB=1, n=5, print_cmd=False, stdout=outfile.replace('.gp', '.log'), stderr=outfile.replace('.gp', '.err'))

#############################################
########## 2. Lift
#############################################

def liftoverJobs():
	for organism in ['human', 'mouse']:
		gp_files = glob.glob('arion/isoseq/s10-liftover.dir/{organism}/gp/*.gp'.format(**locals()))
		genome_str = organism.replace('human', 'hg38').replace('mouse', 'mm10')
		liftover_files = glob.glob('arion/datasets/liftover/{organism}/{genome_str}*.over.chain.gz'.format(**locals()))
		for liftover_file in liftover_files:
			lift_string = os.path.basename(liftover_file).split('.')[0]
			for gp_file in gp_files:
				basename = os.path.basename(gp_file)[:-len('.gp')]
				infiles = [gp_file, liftover_file]
				outfile = 'arion/isoseq/s10-liftover.dir/{organism}/lift/{lift_string}/{basename}-{lift_string}.gp'.format(**locals())
				yield [infiles, outfile]

@follows(convertToGenePred)

@files(liftoverJobs)

def runLiftOver(infiles, outfile):

	# Unmapped
	unmapped_file = outfile.replace('.gp', '_unmapped.gp')

	# Command
	cmd_str = ''' liftOver -genePred {infiles[0]} {infiles[1]} {outfile} {unmapped_file} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['liftover/09-Jul-2019'], W='00:15', GB=1, n=15, print_cmd=False, stdout=outfile.replace('.gp', '.log'), stderr=outfile.replace('.gp', '.err'))

#############################################
########## 3. Merge
#############################################

@follows(runLiftOver)

@collate('arion/isoseq/s10-liftover.dir/*/lift/*/*.gp',
		 regex(r'(.*)/lift/(.*)/(.*)_..(.cds)-(.*.gp)'),
		 r'\1/merged/\2/\3\4-\5')

def mergeLiftOver(infiles, outfile):

	# String
	infiles_str = ' '.join(infiles)

	# Create directory
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Merge
	os.system('cat {infiles_str} > {outfile}'.format(**locals()))

#############################################
########## 4. Filter
#############################################

# @follows(mergeLiftOver)

# @transform('arion/isoseq/s10-liftover.dir/*/merged/*/*.gp',
# 		   regex(r'(.*)/(.*)/(merged)/(.*)(?!d).{1}.gp'),
# 		   add_inputs(r'arion/illumina/s04-alignment.dir/\2/all/gtf/*-all-SJ_filtered.gtf'),
# 		    r'\1/\2/\3/\4_filtered.gp')

# def filterGenePred(infiles, outfile):

# 	# Run
# 	run_r_job('filter_genepred', infiles, outfile, conda_env='env', W='00:05', GB=10, n=1, run_locally=False, stdout=outfile.replace('.gp', '.log'), stderr=outfile.replace('.gp', '.err'))

#############################################
########## 5. Convert
#############################################

@transform(mergeLiftOver,
		   suffix('.gp'),
		   '.gtf')

def convertLiftOver(infile, outfile):

	# Command
	cmd_str = ''' genePredToGtf file {infile} {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['ucsc-utils/2020-03-17'], W='00:15', GB=1, n=15, print_cmd=False, stdout=outfile.replace('.gtf', '.log'), stderr=outfile.replace('.gtf', '.err'))

#############################################
########## 6. Add gene ID
#############################################

@transform(convertLiftOver,
		   regex(r'(.*)/(.*)/(merged)/(.*).gtf'),
		   add_inputs(r'arion/isoseq/s05-talon.dir/\2/gtf/*-SJ_filtered.gtf'),
		#    add_inputs(r'arion/illumina/s04-alignment.dir/\2/all/gtf/*-all-SJ_filtered.gtf'),
		    r'\1/\2/\3/\4-gene_id.gtf')

def addGeneId(infiles, outfile):

	# Run
	run_r_job('add_gene_id', infiles, outfile, conda_env='env', W='00:05', GB=10, n=1, run_locally=False)#, stdout=outfile.replace('.gp', '.log'), stderr=outfile.replace('.gp', '.err'))

#############################################
########## 7. Convert to BigBed
#############################################

def indexJobs():
	gtf_dict = {
		# 'human_ensembl': {
		# 	'gtf': 'arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.102.gtf',
		# 	'chromsize_file': 'arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly_renamed.nochr.chromsizes'
		# },
		# # 'human_all': {
		# # 	'gtf': 'arion/isoseq/s06-cpat.dir/human/gtf/Homo_sapiens.GRCh38.102_talon.cds-all_SJ_filtered.gtf',
		# # 	'chromsize_file': 'arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly_renamed.nochr.chromsizes'
		# # },
		# 'human_novel': {
		# 	'gtf': 'arion/isoseq/s06-cpat.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.cds.novel.gtf',
		# 	'chromsize_file': 'arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly_renamed.nochr.chromsizes'
		# },
		'macaque': {
			'gtf': 'arion/isoseq/s10-liftover.dir/human/merged/hg38ToRheMac10/Homo_sapiens.GRCh38.102_talon-SJ_filtered.cds-hg38ToRheMac10-gene_id.gtf',
			'chromsize_file': 'arion/datasets/reference_genomes/macaque/rheMac10.chrom.sizes'
		},
		'marmoset': {
			'gtf': 'arion/isoseq/s10-liftover.dir/human/merged/hg38ToCalJac4/Homo_sapiens.GRCh38.102_talon-SJ_filtered.cds-hg38ToCalJac4-gene_id.gtf',
			'chromsize_file': 'arion/datasets/reference_genomes/marmoset/calJac4.chrom.sizes'
		},
	}
	for organism, parameter_dict in gtf_dict.items():
		gtf = parameter_dict['gtf']
		infiles = [gtf, parameter_dict['chromsize_file']]
		outfile = os.path.join(os.path.dirname(gtf), 'ucsc_index', os.path.basename(gtf)[:-len('.gtf')]+'.bb')
		yield [infiles, outfile]

@files(indexJobs)

def indexGenomes(infiles, outfile):

	# Parameters
	gff3_file = outfile.replace('.bb', '.gff3')
	genepred_file = outfile.replace('.bb', '.genePred')
	bed_file = outfile.replace('.bb', '.bed')
	input_file = outfile.replace('.bb', '_input.txt')
	index_file = outfile.replace('.bb', '.ix')
	index2_file = outfile.replace('.bb', '.ixx')
	genename_attr = 'gene_name' if 'human' in outfile else 'geneID'

	# Command
	cmd_str = ''' gtfToGenePred -genePredExt {infiles[0]} stdout | sort -k2,2 -k4n,4n > {genepred_file} && \
		genePredToBed {genepred_file} {bed_file} && \
		bedToBigBed -extraIndex=name {bed_file} {infiles[1]} {outfile} && \
		cat {genepred_file} | awk '{{print $1, $12, $1}}' > {input_file} && \
		ixIxx {input_file} {index_file} {index2_file}
	 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['python/2.7.16', 'cufflinks/2.2.1', 'ucsc-utils/2020-03-17'], W='00:15', GB=30, n=1, print_outfile=False, print_cmd=False, stdout=outfile.replace('.bb', '.log'), stderr=outfile.replace('.bb', '.err'))

#######################################################
#######################################################
########## S10. BLAST
#######################################################
#######################################################

#############################################
########## 1. Download genomes
#############################################

def downloadGenomeJobs():
	genomes = pd.read_table('arion/datasets/reference_genomes/ucsc-releases.tsv', comment='#')['genome']
	for genome in genomes:
		yield [None, 'arion/isoseq/s11-blast.dir/genomes/{genome}.fa.gz'.format(**locals())]

@files(downloadGenomeJobs)

def downloadBlastGenomes(infile, outfile):

	# Create directory
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Get genome
	genome = os.path.basename(outfile)[:-len('.fa.gz')]

	# Get URL
	url = 'https://hgdownload.soe.ucsc.edu/goldenPath/{genome}/bigZips/{genome}.fa.gz'.format(**locals())

	# Download
	os.system('cd {outdir} && wget {url}'.format(**locals()))

#############################################
########## 2. Rename genomes
#############################################

@transform(downloadBlastGenomes,
		  suffix('.fa.gz'),
		  '.renamed.fa.gz')

def renameBlastGenomes(infile, outfile):

	# Get genome
	genome = os.path.basename(infile)[:-len('.fa.gz')]

	# Command
	cmd_str = ''' zcat {infile} | sed 's/>/>{genome}_/g' | gzip > {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='01:00', GB=10, n=1, print_outfile=True, print_cmd=False)

#############################################
########## 3. Create database
#############################################

@collate(renameBlastGenomes,
		 regex(r'(.*)/genomes/.*.fa.gz'),
		 r'\1/database/blast_genome_database.db.00.nhr')

def makeBlastDatabase(infiles, outfile):

	# Settings
	basename = outfile.replace('.00.nhr', '')
	merged_fasta = outfile.replace('.db.00.nhr', '.fasta')
	infiles_str = ' '.join(infiles)

	# Command
	cmd_str = ''' zcat {infiles_str} > {merged_fasta} && makeblastdb -in {merged_fasta} -out {basename} -parse_seqids -dbtype nucl && rm {merged_fasta} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, print_cmd=True, modules=['blast/2.9.0+'], W='06:00', GB=30, n=1, stdout=outfile.replace('.db.00.nhr', '.log'), stderr=outfile.replace('.db.00.nhr', '.err'))

#############################################
########## 4. Split
#############################################

# @follows(getFASTA)

# @subdivide('arion/isoseq/s05-talon.dir/*/*_talon.fasta',
		#    regex(r'(.*)/s05-talon.dir/(.*)/(.*).fasta'),
@subdivide('arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.fasta',
		   regex(r'arion/isoseq/s05-talon.dir/(.*)/gtf/(.*).fasta'),
		   r'arion/isoseq/s11-blast.dir/results/\1/fasta/*.fa.*',
		   r'arion/isoseq/s11-blast.dir/results/\1/fasta/\2.fa.1')

def splitBlastFASTA(infile, outfiles, outfileRoot):

	# Get number
	N = 5000 if 'human' in outfileRoot else 100#outfileRoot.split('.')[-1]

	# Get temp filename
	outdir = os.path.dirname(outfileRoot)
	tempfile = os.path.join(outdir, os.path.basename(infile))

	# Command
	cmd_str = ''' cp {infile} {outdir} && gt splitfasta -numfiles {N} {tempfile} && rm {tempfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfileRoot, modules=['genometools/1.5.9'], W='00:15', GB=10, n=1, jobname='splitFASTA_'+outfileRoot.split('/')[-3], stdout=os.path.join(os.path.dirname(outfileRoot), 'job.log'))

#############################################
########## 5. Run BLAST
#############################################

@transform('arion/isoseq/s11-blast.dir/results/human/fasta/*.fasta*',
		   regex(r'(.*)/s11-blast.dir/results/(.*)/fasta/(.*)'),
		   add_inputs(makeBlastDatabase),
		   r'\1/s11-blast.dir/results/\2/split/\3.blast_results.tsv.gz')

def runBLAST(infiles, outfile):

	# Settings
	outfile_tsv = outfile.replace('.gz', '')
	db = infiles[1].replace('.00.nhr', '')

	# Command
	cmd_str = ''' blastn -query {infiles[0]} -db {db} -out {outfile_tsv} -task megablast -outfmt 6 -html -num_threads 1 -evalue 1 && gzip {outfile_tsv}'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, print_outfile=False, modules=['blast/2.9.0+'], W='06:00', GB=10, n=30, q='premium', stdout=outfile.replace('.tsv.gz', '.log'), stderr=outfile.replace('.tsv.gz', '.err'))

#############################################
########## 6. Filter BLAST
#############################################

@transform('arion/isoseq/s11-blast.dir/results/human/split/Homo_sapiens.GRCh38.102_talon-SJ_filtered.fasta.*.blast_results.tsv.gz',
		   regex(r'(.*)/split/(.*).tsv.gz'),
		   r'\1/filtered-75-100/\2_filtered.tsv')

def filterBLAST(infile, outfile):

	# Run
	run_r_job('filter_blast', infile, outfile, conda_env='env', W='00:15', GB=30, n=1, print_outfile=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 6. Filter BLAST
#############################################

@transform('arion/isoseq/s11-blast.dir/results/human/split/Homo_sapiens.GRCh38.102_talon-SJ_filtered.fasta.*.blast_results.tsv.gz',
		   regex(r'(.*)/split/(.*).tsv.gz'),
		   r'\1/filtered-count/\2_filtered.tsv')

def countFilteredBLAST(infile, outfile):

	# Run
	run_r_job('count_filtered_blast', infile, outfile, conda_env='env', W='00:15', GB=30, n=1, print_outfile=True, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#######################################################
#######################################################
########## S. 
#######################################################
#######################################################

#############################################
########## . 
#############################################


##################################################
##################################################
########## Run pipeline
##################################################
##################################################
if __name__ == '__main__':
	cmdline.run(options)
print('Done!')