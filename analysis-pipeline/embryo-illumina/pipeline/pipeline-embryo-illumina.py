#################################################################
#################################################################
############### Embryo Illumina ################
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
import ruffus.cmdline as cmdline
import sys
import os
import json
import glob
import re
# import gtfparse
import pandas as pd
import numpy as np
# from rpy2.robjects import r, pandas2ri
# pandas2ri.activate()

##### 2. LSF #####
# 2.1 Import
sys.path.append('/hpc/users/torred23/pipelines/support')
import lsf

# 2.2 Default parameters
r_source = 'pipeline/scripts/embryo-illumina.R'
py_source = 'pipeline/scripts/EmbryoIllumina.py'
P = 'acc_apollo'
q = 'premium'
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
def run_py_job(func_name, func_input, outfile, W = W, GB = GB, n = n, **kwargs):
	lsf.run_py_job(func_name, func_input, outfile, py_source=py_source, P=P, q=q, W = W, GB = GB, n = n, mkdir=mkdir_val, **kwargs)

##### 3. Custom script imports #####
# 3.1 Python
sys.path.append('pipeline/scripts')
import EmbryoIllumina as S

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
# with open('/sc/arion/projects/GuccioneLab/genome-indices/genome-indices.json') as openfile:
# 	genome_indices = json.load(openfile)

##### 3. Variables #####
# Metadata
organism_metadata = {
	'mouse': 'arion/datasets/qiao/qiao-sample_names.csv',
	'human': 'arion/datasets/human/human-sample_names.csv'
}

# Qiao
# qiao_metadata = 'arion/datasets/qiao/qiao-sample_names.csv'
qiao_illumina = glob.glob('arion/datasets/qiao/rawdata/illumina/*.fastq.gz')

# Human
# human_metadata = 'arion/datasets/human/human-sample_names.csv'
human_illumina = glob.glob('arion/datasets/human/human_embryo_ilmn/*/*/*.fastq.gz')

# All
all_illumina_fastq = glob.glob('arion/illumina/s01-fastq.dir/*/*/*/*.f*q.gz')
trimmed_illumina_fastq = glob.glob('arion/illumina/s01-fastq.dir/*/trimmed/*/*.fq.gz')

# Outlier samples
outlier_sample_file = 'pipeline/outlier_samples.json'
with open(outlier_sample_file) as openfile:
	outlier_samples = json.load(openfile)

# References
reference_dict = {
	'human': {
		'isoseq': {
			'genome_fasta': 'arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa',
			'gtf_filtered': 'arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.gtf',
			'transcript_fasta': 'arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.fasta',
			'cpat_predictions': 'arion/isoseq/s06-cpat.dir/human/Homo_sapiens.GRCh38.102_talon-SJ_filtered-cpat.ORF_prob.best_results.tsv',
			'pfam_predictions': 'arion/isoseq/s07-pfam.dir/human/human-translated_pfam.tsv',
			'gtf_cds': 'arion/isoseq/s06-cpat.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.cds.gtf'
		}
	},
	'mouse': {
		'isoseq': {
			'genome_fasta': 'arion/datasets/reference_genomes/mouse/Mus_musculus.GRCm38.dna_sm.primary_assembly.fa',
			'gtf_filtered': 'arion/isoseq/s05-talon.dir/mouse/gtf/Mus_musculus.GRCm38.102_talon-SJ_filtered.gtf',
		}
	}
}

# Comparisons
comparison_file = 'pipeline/comparisons.json'
with open(comparison_file) as openfile:
	comparison_dict = json.load(openfile)

#######################################################
#######################################################
########## S1. FASTQ
#######################################################
#######################################################

#############################################
########## 1. Link
#############################################

def linkJobs():

	# Read metadata
	mouse_dataframe = pd.read_csv(organism_metadata['mouse']).query('Platform == "ILLUMINA"')
	human_dataframe = pd.read_csv(organism_metadata['human']).query('Platform == "ILLUMINA"')
	
	# Samples
	sample_dict = {
		'mouse': [{'sample_name': rowData['sample_name'], 'fastq': [x for x in qiao_illumina if rowData['Run'] in x]} for  index, rowData in mouse_dataframe.iterrows()],
		'human': [{'sample_name': rowData['sample_name'], 'fastq': [x for x in human_illumina if x.split('/')[-3] == rowData['sample_id']]} for index, rowData in human_dataframe.iterrows()]
	}

	# Loop
	for organism, samples in sample_dict.items():
		for sample in samples:
			for fastq_file in sample['fastq']:

				# Get FASTQ name
				fastq_basename = os.path.basename(fastq_file)

				# Add infile path
				infile = os.path.join(os.getcwd(), fastq_file)

				# Get outfile
				outfile = 'arion/illumina/s01-fastq.dir/{organism}/raw/{sample_name}/{fastq_basename}'.format(**locals(), **sample)

				# Yield
				yield [infile, outfile]

@files(linkJobs)

def linkFASTQ(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)

	# Create
	os.system('mkdir -p {outdir} && ln -s {infile} {outfile}'.format(**locals()))

#############################################
########## 2. Adapter trimming
#############################################

def trimJobs():
	for sample_path in glob.glob('arion/illumina/s01-fastq.dir/*/raw/*'):
		infiles = glob.glob(os.path.join(sample_path, '*'))
		infiles.sort()
		outdir = sample_path.replace('/raw/', '/trimmed/')
		yield [infiles, outdir]

@follows(linkFASTQ)

@files(trimJobs)

def trimIlluminaAdapters(infiles, outdir):

	# Command
	cmd_str = '''trim_galore --nextera --paired --cores 6 --output_dir {outdir} {infiles[0]} {infiles[1]}'''.format(**locals())

	# Run
	run_job(cmd_str, outdir, modules=['trim_galore/0.6.6'], W='06:00', GB=6, n=6, stdout=os.path.join(outdir, 'job.log'), stderr=os.path.join(outdir, 'job.err'))

#######################################################
#######################################################
########## S2. QC
#######################################################
#######################################################

#############################################
########## 1. FASTQC
#############################################

@follows(linkFASTQ, trimIlluminaAdapters)

@transform(all_illumina_fastq,
		   regex(r'(.*)/s01-fastq.dir/(.*).f.*q.gz'),
		   r'\1/s02-fastqc.dir/\2_fastqc.html')

def runFastQC(infile, outfile):

	# Command
	cmd_str = '''fastqc --outdir=$(dirname {outfile}) {infile}'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['fastqc/0.11.8'], W='03:00', GB=12, n=1, print_outfile=False)

#############################################
########## 2. MultiQC
#############################################

# @follows(runFastQC)

def qcJobs():
	filelist = [
		['arion/illumina/s02-fastqc.dir/human/raw', 'arion/illumina/multiqc/human_fastqc/multiqc_report.html'],
		['arion/illumina/s02-fastqc.dir/mouse/raw', 'arion/illumina/multiqc/mouse_fastqc/multiqc_report.html'],
		['arion/illumina/s02-fastqc.dir/human/trimmed', 'arion/illumina/multiqc/human_fastqc_trimmed/multiqc_report.html'],
		['arion/illumina/s02-fastqc.dir/mouse/trimmed', 'arion/illumina/multiqc/mouse_fastqc_trimmed/multiqc_report.html'],
		['arion/illumina/s03-alignment.dir/human', 'arion/illumina/multiqc/human_alignment_trimmed/multiqc_report.html'],
		['arion/illumina/s03-alignment.dir/mouse', 'arion/illumina/multiqc/mouse_alignment_trimmed/multiqc_report.html'],
		['arion/illumina/s03-alignment.dir-old/human/isoseq/all', 'arion/illumina/multiqc/human_test_qc/multiqc_report.html'],
	]
	for files in filelist:
		yield [files[0], files[1]]

@files(qcJobs)

def runMultiQC(infile, outfile):

	# Command
	cmd_str = 'multiqc --outdir $(dirname {outfile}) {infile}'.format(**locals())

	# Run
	if not os.path.exists(outfile):
		run_job(cmd_str, outfile, conda_env='env', W="01:00", GB=10, n=1, print_outfile=False, run_locally=False, stdout=outfile.replace('.html', '.log'))

#######################################################
#######################################################
########## S3. Alignment and quantification
#######################################################
#######################################################

#############################################
########## 1. STAR index
#############################################

def starIndexJobs():
	for organism, organism_references in reference_dict.items():
		for source, reference_files in organism_references.items():
			outfile = 'arion/illumina/s03-alignment.dir/{organism}/{source}/STAR/index'.format(**locals())
			yield [reference_files, outfile]

# @follows(runFastQC)

@files(starIndexJobs)

def buildStarIndex(infiles, outfile):

	# Command
	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {genome_fasta} --sjdbGTFfile {gtf_filtered} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals(), **infiles)

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='02:00', GB=5, n=15, ow=True, print_cmd=True, stdout=os.path.join(outfile, 'job.log'), jobname='_'.join(outfile.split('/')[-4:]), wait=False)

#############################################
########## 2. STAR junctions
#############################################

def starJunctionJobs():
	fastq_dataframe = pd.DataFrame([{'fastq': x, 'organism': x.split('/')[-4], 'sample_name': x.split('/')[-2]} for x in trimmed_illumina_fastq]).sort_values('fastq').groupby(['organism', 'sample_name'])['fastq'].apply(list).reset_index()
	for organism, sample_dataframe in fastq_dataframe.groupby('organism'):
		fastq_dict = sample_dataframe.drop('organism', axis=1).set_index('sample_name')['fastq'].to_dict()
		for sample_name, fastq_files in fastq_dict.items():
			if sample_name not in outlier_samples[organism]:
				for source in ['isoseq']:
					star_index = 'arion/illumina/s03-alignment.dir/{organism}/{source}/STAR/index'.format(**locals())
					outfile = 'arion/illumina/s03-alignment.dir/{organism}/{source}/STAR/pass1/{sample_name}/{sample_name}-SJ.out.tab'.format(**locals())
					yield [(fastq_files, star_index), outfile]

# @follows(linkFASTQ, buildStarIndex)

@files(starJunctionJobs)

def getStarJunctions(infiles, outfile):

	# Split
	fastq_files, star_index = infiles

	# Prefix
	prefix = outfile[:-len('SJ.out.tab')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_files[0]} {fastq_files[1]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 100 \
		--outSAMtype None'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, print_outfile=False, W="06:00", GB=5, n=15, modules=['star/2.7.5b'], stdout=outfile.replace('-SJ.out.tab', '_job.log'))

#############################################
########## 3. STAR BAM
#############################################

def starJobs():
	fastq_dataframe = pd.DataFrame([{'fastq': x, 'organism': x.split('/')[-4], 'sample_name': x.split('/')[-2]} for x in trimmed_illumina_fastq]).sort_values('fastq').groupby(['organism', 'sample_name'])['fastq'].apply(list).reset_index()
	for organism, sample_dataframe in fastq_dataframe.groupby('organism'):
		fastq_dict = sample_dataframe.drop('organism', axis=1).set_index('sample_name')['fastq'].to_dict()
		for sample_name, fastq_files in fastq_dict.items():
			if sample_name not in outlier_samples[organism]:
				for source in ['isoseq']:
					star_index = 'arion/illumina/s03-alignment.dir/{organism}/{source}/STAR/index'.format(**locals())
					sj_files = glob.glob('arion/illumina/s03-alignment.dir/{organism}/{source}/STAR/pass1/*/*-SJ.out.tab'.format(**locals()))
					outfile = 'arion/illumina/s03-alignment.dir/{organism}/{source}/STAR/pass2/{sample_name}/{sample_name}-Aligned.sortedByCoord.out.bam'.format(**locals())
					yield [(fastq_files, star_index, sj_files), outfile]

@follows(getStarJunctions)

@files(starJobs)

def runStar(infiles, outfile):

	# Split
	fastq_files, star_index, sj_files = infiles

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
	sj_files_str = ' '.join(sj_files)

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_files[0]} {fastq_files[1]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 5000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="06:00", GB=10, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.log'))

#############################################
########## 4. RSEM index
#############################################

def rsemIndexJobs():
	for organism, organism_references in reference_dict.items():
		for source, reference_files in organism_references.items():
			gtf_name = os.path.basename(reference_files['gtf_filtered'])[:-len('.gtf')]
			outfile = 'arion/illumina/s03-alignment.dir/{organism}/{source}/RSEM/index/{gtf_name}.idx.fa'.format(**locals())
			yield [reference_files, outfile]

@files(rsemIndexJobs)

def createRsemIndex(infiles, outfile):

	# Command
	basename = outfile[:-len('.idx.fa')]
	cmd_str = ''' rsem-prepare-reference --gtf {gtf_filtered} --num-threads 10 {genome_fasta} {basename} '''.format(**locals(), **infiles)

	# Run
	run_job(cmd_str, outfile, W="00:30", GB=5, n=3, modules=['rsem/1.3.3'], print_cmd=True, stdout=basename+'.log', stderr=basename+'.err')

#############################################
########## 5. RSEM expression
#############################################

# @follows(runStar, createRsemIndex)

@transform('arion/illumina/s03-alignment.dir/*/*/STAR/pass2/*/*-Aligned.toTranscriptome.out.bam',
		   regex(r'(.*)/STAR/.*/(.*)-Aligned.toTranscriptome.out.bam'),
		   add_inputs(r'\1/RSEM/index/*.idx.fa'),
		   r'\1/RSEM/results/\2/\2.isoforms.results')

def runRsem(infiles, outfile):

	# Variables
	prefix = outfile[:-len('.isoforms.results')]
	reference_name = infiles[1][:-len('.idx.fa')]

	# Command
	cmd_str = '''rsem-calculate-expression \
		--alignments \
		--strandedness none \
		--paired-end \
		--estimate-rspd \
		--num-threads 200 \
		{infiles[0]} \
		{reference_name} \
		{prefix} > {prefix}.rsem.log && \
		rsem-plot-model {prefix} {prefix}.quant.pdf '''.format(**locals())
		# --calc-ci \

	# Run
	run_job(cmd_str, outfile, W="06:00", GB=3, n=25, modules=['rsem/1.3.3'], print_cmd=True, stdout=outfile.replace('.isoforms.results', '.log'), stderr=outfile.replace('.isoforms.results', '.err'))

#############################################
########## 6. Filter BAM
#############################################

@transform('arion/illumina/s03-alignment.dir/human/isoseq/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
		   suffix('.bam'),
		   '.filtered.bam')

def filterBam(infile, outfile):

	# Files
	outname = outfile.replace('.bam', '')

	# Command
	cmd_str = ''' sambamba view --with-header --nthreads 30 --format bam --filter "not unmapped and not duplicate and not secondary_alignment and mapping_quality >= 30" {infile} > {outfile} && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='03:00', n=5, GB=5, modules=['samtools/1.9', 'sambamba/0.5.6'], print_outfile=False)

#############################################
########## 4. Junction counts
#############################################

# def sjCountJobs():
# 	for organism, organism_references in reference_dict.items():
# 		for source, reference_files in organism_references.items():
# 			infiles = [reference_files['talon_junctions']] + glob.glob('arion/illumina/s03-junctions.dir/{organism}/isoseq/STAR/pass2/*/*-SJ.out.tab'.format(**locals()))
# 			outfile = 'arion/illumina/s03-junctions.dir/{organism}/isoseq/STAR/{organism}-{source}-junction_counts.tsv'.format(**locals())
# 			yield [infiles, outfile]

# @follows(runStar)

# @files(sjCountJobs)

# def getJunctionCounts(infiles, outfile):

# 	# Run
# 	run_r_job('get_junction_counts', infiles, outfile, run_locally=False, conda_env='env', W='02:00', GB=50, n=1, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

# #######################################################
# #######################################################
# ########## S4. Alignment
# #######################################################
# #######################################################

# #############################################
# ########## 1. Filter GTF
# #############################################

# def filterJobs():
# 	for organism, comparisons in comparison_dict.items():
# 		source = 'isoseq'
# 		reference_files = reference_dict[organism][source]
# 		gtf = reference_files['gtf']
# 		talon_abundance_file = reference_files['talon_abundance']
# 		sj_counts_file = 'arion/illumina/s03-junctions.dir/{organism}/{source}/STAR/{organism}-{source}-junction_counts.tsv'.format(**locals())
# 		infiles = [gtf, talon_abundance_file, sj_counts_file, outlier_sample_file]
# 		for comparison in comparisons+['all']:
# 			if comparison == 'all':
# 				gtf_prefix = os.path.basename(gtf)[:-len('.gtf')]
# 				comparison_string = '_vs_'.join(comparison) if comparison != 'all' else comparison
# 				outfile = 'arion/illumina/s03-alignment.dir/{organism}/{comparison_string}/gtf/{gtf_prefix}-{comparison_string}-SJ_filtered.gtf'.format(**locals())
# 				yield [infiles, outfile, comparison]

# # @follows(getJunctionCounts)

# @files(filterJobs)

# def filterGTF(infiles, outfile, comparison):

# 	# Run
# 	run_r_job('filter_gtf', infiles, outfile, additional_params=comparison, run_locally=False, conda_env='env', W='00:10', GB=10, n=1, stdout=outfile.replace('.gtf', '.log'), stderr=outfile.replace('.gtf', '.err'))

# #############################################
# ########## 2. STAR index
# #############################################

# @transform(filterGTF,
# 		   regex(r'(.*.dir)/(.*?)/(.*)/gtf/.*.gtf'),
# 		   add_inputs(r'arion/datasets/reference_genomes/\2/*.dna_sm.primary_assembly.fa'),
# 		   r'\1/\2/\3/STAR/index')

# def buildStarIndexFiltered(infiles, outfile):

# 	# Command
# 	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {infiles[1]} --sjdbGTFfile {infiles[0]} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='02:00', GB=3, n=15, ow=True, print_cmd=False, stdout=os.path.join(outfile, 'job.log'), jobname='_'.join(outfile.split('/')[-4:]), wait=False)

# # find arion/illumina/s03-alignment.dir/*/isoseq/*/STAR/index -name "job.log" | jsc

# #############################################
# ########## 3. RSEM index
# #############################################

# # @transform(filterGTF,
# # 		   regex(r'(arion)/(.*.dir)/(.*?)/(.*)/gtf/(.*).gtf'),
# # 		   add_inputs(r'\1/datasets/reference_genomes/\3/*.dna_sm.primary_assembly.fa'),
# # 		   r'\1/\2/\3/\4/RSEM/index/\5.idx.fa')

# # def createRsemIndex(infiles, outfile):

# # 	# Command
# # 	basename = outfile[:-len('.idx.fa')]
# # 	cmd_str = ''' rsem-prepare-reference --gtf {infiles[0]} --num-threads 10 {infiles[1]} {basename} '''.format(**locals())

# # 	# Run
# # 	run_job(cmd_str, outfile, W="00:30", GB=5, n=3, modules=['rsem/1.3.3'], print_cmd=False, stdout=basename+'.log', stderr=basename+'.err')

# # find arion/illumina/s03-alignment.dir/*/*/*/RSEM/index -name "*.log" | jsc

# #############################################
# ########## 4. STAR
# #############################################

# def starFilteredJobs():

# 	# Read samples
# 	fastq_dataframe = pd.DataFrame([{'fastq': x, 'organism': x.split('/')[-4], 'sample_name': x.split('/')[-2]} for x in trimmed_illumina_fastq]).sort_values('fastq').groupby(['organism', 'sample_name'])['fastq'].apply(list).reset_index()

# 	# Loop through organisms
# 	for organism, sample_dataframe in fastq_dataframe.groupby('organism'):
# 		fastq_dict = sample_dataframe.drop('organism', axis=1).set_index('sample_name')['fastq'].to_dict()

# 		# Loop through comparisons
# 		for comparison in comparison_dict[organism]+['all']:
# 			comparison_string = '_vs_'.join(comparison) if comparison != 'all' else comparison

# 			# Loop throush samples
# 			for sample_name, fastq_files in fastq_dict.items():
# 				cell_type = sample_name.split('_')[1].replace('2PN', '1C')

# 				# Select samples for each comparison
# 				# if cell_type==comparison[0] or cell_type==comparison[1] or comparison == 'all':
# 				if comparison == 'all' and sample_name not in outlier_samples[organism]:
# 					star_index = 'arion/illumina/s03-alignment.dir/{organism}/{comparison_string}/STAR/index'.format(**locals())
# 					sj_files = glob.glob('arion/illumina/s03-junctions.dir/{organism}/isoseq/STAR/pass1/*/*-SJ.out.tab'.format(**locals()))
# 					outfile = 'arion/illumina/s03-alignment.dir/{organism}/{comparison_string}/STAR/pass2/{sample_name}/{sample_name}-Aligned.sortedByCoord.out.bam'.format(**locals())
# 					yield [(fastq_files, star_index, sj_files), outfile]

# @follows(getStarJunctions, buildStarIndexFiltered)

# @files(starFilteredJobs)

# def runStarFiltered(infiles, outfile):

# 	# Split
# 	fastq_files, star_index, sj_files = infiles

# 	# Variables
# 	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
# 	sj_files_str = ' '.join(sj_files)

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {star_index} \
# 		--readFilesIn {fastq_files[0]} {fastq_files[1]} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 32 \
# 		--sjdbFileChrStartEnd {sj_files_str} \
# 		--limitSjdbInsertNsj 5000000 \
# 		--quantMode TranscriptomeSAM GeneCounts \
# 		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="06:00", GB=10, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.log'))

# # find arion/illumina/s03-alignment.dir/human/*/STAR/pass2 -name "*.log" | jsc

# #############################################
# ########## 5. RSEM expression
# #############################################

# @follows(runStarFiltered, createRsemIndex)

# # @transform('arion/illumina/s03-alignment.dir/human/2C_vs_4C/STAR/pass2/human_4C_B3_9/human_4C_B3_9-Aligned.toTranscriptome.out.bam',
# @transform('arion/illumina/s03-alignment.dir/*/*/STAR/pass2/*/*-Aligned.toTranscriptome.out.bam',
# 		   regex(r'(.*)/STAR/.*/(.*)-Aligned.toTranscriptome.out.bam'),
# 		   add_inputs(r'\1/RSEM/index/*.idx.fa'),
# 		   r'\1/RSEM/results/\2/\2.isoforms.results')

# def runRsem2(infiles, outfile):

# 	# Variables
# 	prefix = outfile[:-len('.isoforms.results')]
# 	reference_name = infiles[1][:-len('.idx.fa')]

# 	# Command
# 	cmd_str = '''rsem-calculate-expression \
# 		--alignments \
# 		--strandedness none \
# 		--paired-end \
# 		--estimate-rspd \
# 		--num-threads 200 \
# 		{infiles[0]} \
# 		{reference_name} \
# 		{prefix} > {prefix}.rsem.log && \
# 		rsem-plot-model {prefix} {prefix}.quant.pdf '''.format(**locals())
# 		# --calc-ci \

# 	# Run
# 	run_job(cmd_str, outfile, W="06:00", GB=2, n=25, modules=['rsem/1.3.3'], print_cmd=True, stdout=outfile.replace('.isoforms.results', '.log'), stderr=outfile.replace('.isoforms.results', '.err'))

#############################################
########## 6. Create BigWig
#############################################

# @transform(runStarFiltered,
# 		   suffix('-Aligned.sortedByCoord.out.bam'),
# 		   '.bw')

# def createBigWig(infile, outfile):

# 	# Command
# 	cmd_str = """bamCoverage --outFileFormat=bigwig --skipNonCoveredRegions --numberOfProcessors=50 --normalizeUsing RPKM -b {infile} -o {outfile}""".format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', W='05:00', n=7, GB=6)

#############################################
########## 7. Filtered FASTA
#############################################

# @transform('arion/illumina/s03-alignment.dir/*/all/gtf/*SJ_filtered.gtf',
# 		   regex(r'(.*.dir)/(.*?)/(.*)/(.*).gtf'),
# 		   add_inputs(r'arion/datasets/reference_genomes/\2/*dna_sm.primary_assembly.fa'),
# 		   r'\1/\2/\3/SQANTI3/\4.fasta')

# def getFilteredFasta(infiles, outfile):

# 	# Command
# 	cmd_str = ''' gffread {infiles[0]} -g {infiles[1]} -w {outfile} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['gff/2021-02'], W='00:30', GB=10, n=1, stdout=outfile.replace('.fasta', '_fasta.log'), stderr=outfile.replace('.fasta', '_fasta.err'))

#############################################
########## 8. Run SQANTI
#############################################

# @transform('arion/illumina/s03-alignment.dir/human/all/gtf/Homo_sapiens.GRCh38.102_talon-all-SJ_filtered.gtf',
# 		   regex(r'(.*.dir)/(.*?)/(.*)/(.*).gtf'),
# 		   add_inputs(r'arion/datasets/reference_genomes/\2/*.gtf', r'arion/datasets/reference_genomes/\2/*.dna_sm.primary_assembly.fa', '/sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.1/data/polyA_motifs/mouse_and_human.polyA_motif.txt', [r'arion/illumina/s03-alignment.dir/\2/all/STAR/pass2/*/*-SJ.out.tab']),
# 		   r'\1/\2/\3/SQANTI3_test5/\4-SQANTI3_report.pdf')

# def runSqanti(infiles, outfile):

# 	# Prepare
# 	sj_files = ','.join(infiles[4])
# 	dirname = os.path.dirname(outfile)
# 	basename = os.path.basename(outfile)[:-len('_report.pdf')]

# 	# Command
# 			# --polyA_motif_list {infiles[3]} \
# 			# --cpus 10 \
# 	cmd_str = ''' export PYTHONPATH=$PYTHONPATH:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.1/cDNA_Cupcake/ && export PYTHONPATH=$PYTHONPATH:/sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.1/cDNA_Cupcake/sequence/ && \
# 		python /sc/arion/work/torred23/libraries/SQANTI/SQANTI3-4.1/sqanti3_qc.py {infiles[0]} {infiles[1]} {infiles[2]} \
# 			--skipORF \
# 			-c {sj_files} \
# 			--report both \
# 			-d {dirname} \
# 			-o {basename}
# 	'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='SQANTI3_v4.1', W='02:00', n=3, GB=10, print_cmd=True)

# find arion/illumina/s03-alignment.dir/*/all/gtf/SQANTI3 -name "SQANTI3" | xargs rm -r

#######################################################
#######################################################
########## S4. Expression
#######################################################
#######################################################

#############################################
########## 1. Prepare metadata
#############################################

@collate(runRsem,
		 regex(r'(.*)/s03-alignment.dir/(.*)/(.*)/RSEM/.*/.*.isoforms.results'),
		 r'\1/s04-expression.dir/\2/\2-sample_metadata.txt')

def prepareSampleMetadata(infiles, outfile):
	
	# Create directory
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Prepare dict
	sample_dict = [{'files': x, 'names': x.split('/')[-2]} for x in infiles]

	# Convert to dataframe
	sample_dataframe = pd.DataFrame(sample_dict)

	# Add more columns
	organism = outfile.split('/')[-2]
	sample_dataframe['organism'] = [organism for x in sample_dataframe['names']]
	sample_dataframe['cell_type'] = [x.split('_')[1].replace('2PN', '1C') for x in sample_dataframe['names']]

	# Read metadata
	metadata_dataframe = pd.read_csv(organism_metadata[organism]).query('Platform == "ILLUMINA"').rename(columns={'sample_name': 'names'}).drop(['Platform', 'count'], axis=1)

	# Fix organism specific columns
	if organism == 'human':
		metadata_dataframe = metadata_dataframe.drop(['Group', 'batch'], axis=1)
		sample_dataframe['batch'] = [x.split('_')[2] for x in sample_dataframe['names']]
	elif organism == 'mouse':
		metadata_dataframe['strain'] = [x.replace(' ', '').replace('/', '_').replace('x', '_') for x in metadata_dataframe['Strain']]
		metadata_dataframe = metadata_dataframe.drop(['Run', 'Cell_type', 'sample_id', 'Strain'], axis=1)
		sample_dataframe['replicate'] = [x.split('_')[-1] for x in sample_dataframe['names']]

	# Merge
	merged_dataframe = sample_dataframe.merge(metadata_dataframe, on='names')

	# Write
	merged_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Aggregate
#############################################

@transform(prepareSampleMetadata,
		   regex(r'(.*)/(s04-expression.dir)/(.*)/.*-sample_metadata.txt'),
		   add_inputs(r'arion/isoseq/s05-talon.dir/\3/gtf/*-SJ_filtered.gtf'),
		   r'\1/\2/\3/\3-counts.rda')

def aggregateCounts(infiles, outfile):

	# Run
	run_r_job('aggregate_counts', infiles, outfile, conda_env='env', run_locally=False, W='00:30', GB=30, n=1, wait=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

# find arion/illumina/s04-expression.dir/*/all -name "*counts*" | xargs rm

#############################################
########## 3. Transcript TPM
#############################################

@transform(aggregateCounts,
		   suffix('counts.rda'),
		   'transcript_tpm.txt')

def getTranscriptTPM(infile, outfile):

	# Run
	run_r_job('get_transcript_tpm', infile, outfile, conda_env='env', run_locally=False, ow=False)

#############################################
########## 4. Gene counts
#############################################

@transform(aggregateCounts,
		   suffix('counts.rda'),
		   'gene_normalized_counts.tsv')

def getGeneExpression(infile, outfile):

	# Run
	run_r_job('get_gene_expression', infile, outfile, conda_env='env', run_locally=False, ow=False)

#############################################
########## 5. Get size factors
#############################################

@transform(aggregateCounts,
		   suffix('-counts.rda'),
		   '-size_factors.tsv')

def getSizeFactors(infile, outfile):

	# Run
	run_r_job('get_size_factors', infile, outfile, conda_env='env', modules=[], W='00:15', GB=20, n=1, run_locally=False, print_outfile=False, print_cmd=False)

	# # Read size factor
	# normalization_dict = pd.read_table(infiles[1], index_col='sample_name')['size_factor_reciprocal'].to_dict()
	# size_factor = normalization_dict[os.path.basename(outfile).split('-')[0]]

	# # Command
	# cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=3 --skipNonCoveredRegions --numberOfProcessors=48 --scaleFactor {size_factor} -b {infiles[0]} -o {outfile}""".format(**locals())

	# # Run
	# run_job(cmd_str, outfile, conda_env='env', W='03:00', n=8, GB=4, print_cmd=False)

#############################################
########## 6. Create BigWig
#############################################

@transform('arion/illumina/s03-alignment.dir/human/isoseq/STAR/pass2/*/*-Aligned.sortedByCoord.out.filtered.bam',
		   suffix('.bam'),
		   '.bw')

def createBigWig(infile, outfile):

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=3 --normalizeUsing CPM --skipNonCoveredRegions --numberOfProcessors=48 -b {infile} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='03:00', n=8, GB=4, print_cmd=False)

#############################################
########## 6. Create scaled BigWig
#############################################

@transform('arion/illumina/s03-alignment.dir/human/isoseq/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
		   regex(r'(.*.dir)/(.*?)/(.*)-Aligned.sortedByCoord.out.bam'),
		   add_inputs(r'arion/illumina/s04-expression.dir/\2/\2-size_factors.tsv'),
		   r'\1/\2/\3-scaled_3.bw')

def createScaledBigWig(infiles, outfile):

	# Read size factor
	normalization_dict = pd.read_table(infiles[1], index_col='sample_name')['size_factor_reciprocal'].to_dict()
	size_factor = normalization_dict[os.path.basename(outfile).split('-')[0]]

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=3 --skipNonCoveredRegions --numberOfProcessors=48 --scaleFactor {size_factor} -b {infiles[0]} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='03:00', n=8, GB=4, print_cmd=False)

#############################################
########## 7. Merge scaled BigWig
#############################################

# @collate(createScaledBigWig,
@collate('arion/illumina/s03-alignment.dir/human/isoseq/STAR/pass2/*/*-scaled_3.bw',
		#  regex(r'(.*)/s03-alignment.dir/(.*)/all/STAR/pass2/.*?_(morula)_.*/.*-scaled_3.bw'),
		 regex(r'(.*)/s03-alignment.dir/(.*)/isoseq/STAR/pass2/.*?_(.*?)_.*/.*.bw'),
		 add_inputs('arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.chrom.sizes'),
		 r'\1/s04-expression.dir/\2/scaled_bw/\2-\3.bw')

def mergeScaledBigWig(infiles, outfile):
	
	# Files
	wig_file = outfile.replace('.bw', '.wig')
	bedgraph_file = outfile.replace('.bw', '.bedgraph')
	infiles_str = ' '.join([x[0] for x in infiles])

	# Command
	cmd_str = """ wiggletools mean {infiles_str} | sed '/^KI.*/{{s///;q;}}' > {wig_file} && wigToBigWig {wig_file} {infiles[0][1]} {outfile} && rm {wig_file} """.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['wiggletools/1.2', 'ucsc-utils/2020-03-17'], W='00:30', n=1, GB=10, print_cmd=True, stdout=outfile.replace('.bw', '.log'), stderr=outfile.replace('.bw', '.err'))

# wiggletools write_bg - arion/illumina/s04-expression.dir/human/all/scaled_bw/human-morula.wig | bedSort | head

# find arion/illumina/s03-alignment.dir/human/all/STAR/pass2 -name "*scaled.bw"

#############################################
########## 7. Merge scaled BigWig
#############################################

# @collate(createScaledBigWig,
@collate('arion/illumina/s03-alignment.dir/human/isoseq/STAR/pass2/*/*.filtered.bw',
		 regex(r'(.*)/s03-alignment.dir/(.*)/isoseq/STAR/pass2/.*?_(.*?)_.*/.*.bw'),
		 add_inputs('arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.chrom.sizes'),
		 r'\1/s04-expression.dir/\2/filtered_bw/\2-\3.bw')

def mergeBigWig(infiles, outfile):
	
	# Files
	wig_file = outfile.replace('.bw', '.wig')
	bedgraph_file = outfile.replace('.bw', '.bedgraph')
	infiles_str = ' '.join([x[0] for x in infiles])

	# Command
	cmd_str = """ wiggletools mean {infiles_str} | sed '/^KI.*/{{s///;q;}}' > {wig_file} && wigToBigWig {wig_file} {infiles[0][1]} {outfile} && rm {wig_file} """.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['wiggletools/1.2', 'ucsc-utils/2020-03-17'], W='00:30', n=1, GB=10, print_cmd=False, stdout=outfile.replace('.bw', '.log'), stderr=outfile.replace('.bw', '.err'))

# wiggletools write_bg - arion/illumina/s04-expression.dir/human/all/scaled_bw/human-morula.wig | bedSort | head

# find arion/illumina/s03-alignment.dir/human/all/STAR/pass2 -name "*scaled.bw"

#############################################
########## 8. Cluster novel genes
#############################################

@transform('arion/illumina/s04-expression.dir/human/human-counts.rda',
		   suffix('counts.rda'),
		   'novel_gene_clusters.rda')

def clusterNovelGenes(infile, outfile):

	# Run
	run_r_job('cluster_novel_genes', infile, outfile, modules=['R/4.0.3'], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 8. Run VIPER
#############################################

@transform('arion/illumina/s04-expression.dir/human/human-counts.rda',
		   regex(r'(.*)/(.*)/.*counts.rda'),
		   add_inputs(r'arion/isoseq/s05-talon.dir/\2/gtf/*-SJ_filtered.gtf'),
		   r'\1/\2/\2-tf_activity.rda')

def runViper(infiles, outfile):

	# Run
	run_r_job('get_tf_activity', infiles, outfile, modules=['R/4.0.3'], W='00:15', GB=15, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 8. Cluster RBPs
#############################################

@transform('arion/illumina/s04-expression.dir/human/human-gene_normalized_counts.tsv',
		   regex(r'(.*).tsv'),
		   add_inputs('arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.102.gtf'),
		   r'\1_splicing_genes_mfuzz.rda')

def clusterRbps(infiles, outfile):

	# Run
	run_r_job('cluster_rbps', infiles, outfile, modules=['R/4.0.3'], W='00:15', GB=15, n=1, run_locally=True)#, print_outfile=False, print_cmd=False)

#######################################################
#######################################################
########## S5. Differential expression
#######################################################
#######################################################

#############################################
########## 1. DESeq2
#############################################

# @follows(aggregateCounts)

@subdivide('arion/illumina/s04-expression.dir/*/*-counts.rda',
		   regex(r'(.*)/s04-expression.dir/(.*)/.*-counts.rda'),
		   add_inputs(comparison_file),
		   r'\1/s05-differential_expression.dir/\2/\2-*-deseq.tsv',
		   r'\1/s05-differential_expression.dir/\2/\2-{comparison[1]}_vs_{comparison[2]}-{feature}-deseq_no_outlier_cutoff.tsv')

def runDESeq2(infiles, outfiles, outfileRoot):

	# Loop
	for feature in ['gene', 'transcript']:

		# Get info
		logname = os.path.join(os.path.dirname(outfileRoot), 'deseq-'+feature+'.log')
		jobname = '_'.join(logname[:-len('.log')].split('/')[-3:])
		
		# Run
		run_r_job('run_deseq2', infiles, outfileRoot, conda_env='env', W='24:00', GB=25, n=1, additional_params=feature, ow=False, stdout=logname, jobname=jobname)

# #############################################
# ########## 2. Excel
# #############################################

# @follows(runDESeq2)

@collate('arion/illumina/s05-differential_expression.dir/human/human-*-gene-deseq.tsv',
		 regex(r'(.*human)-.*.tsv'),
		 add_inputs('arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.gtf'),
		 'arion/figures/supplementary_tables/SupplementaryTable7.xlsx')
		#  r'\1-merged.xlsx')

def createExcel(infiles, outfile):

	# Split infiles
	gtf_file = infiles[0][1]
	infiles = [x[0] for x in infiles]

	# Read differential expression
	dataframes = {x.split('-')[-3]: pd.read_table(x).drop(['lfcSE', 'stat'], axis=1) for x in infiles}

	# Read gene names
	gene_dataframe = gtfparse.read_gtf(gtf_file)[['gene_id', 'gene_name']].drop_duplicates()

	# Read Dalit's genelist
	# dalit_genelist = pd.read_csv('arion/datasets/dalit/Developmental_genelist_by_Dalit.csv').set_index('Gene_Symbol')['Classification'].to_dict()

	# Initialize ExcelWriter
	with pd.ExcelWriter(outfile) as writer:

		# Loop
		for comparison in ['1C_vs_2C', '2C_vs_4C', '4C_vs_8C', '8C_vs_morula', 'morula_vs_blastocyst']:

			# Get data
			dataframe = gene_dataframe.merge(dataframes[comparison], on='gene_id', how='right').sort_values('padj')

			# Split
			comparison_groups = comparison.split('_vs_')

			# Significance
			# dataframe['gene_biotype'] = [rowData['gene_biotype'].replace('_', ' ') if rowData['gene_biotype'] != 'not_available' else rowData['gene_category'].replace('_', ' ') for index, rowData in dataframe.iterrows()]

			# Significance
			dataframe['significant'] = [rowData['padj'] < 0.05 and abs(rowData['log2FoldChange']) > 1 for index, rowData in dataframe.iterrows()]
			dataframe['differential_expression'] = ['Not significant' if not rowData['significant'] else 'Up in '+comparison_groups[1] if rowData['log2FoldChange'] > 0 else 'Down in '+comparison_groups[1] for index, rowData in dataframe.iterrows()]

			# Round
			dataframe['baseMean'] = dataframe['baseMean'].round(1)
			dataframe['log2FoldChange'] = dataframe['log2FoldChange'].round(1)
			dataframe['pvalue'] = ['{:.2e}'.format(x) for x in dataframe['pvalue']]
			dataframe['padj'] = ['{:.2e}'.format(x) for x in dataframe['padj']]

			# Dalit geneset
			# dataframe['Developmental category'] = [dalit_genelist.get(x, '').replace('_', ' ').title() for x in dataframe['gene_name']]
			
			# Rename
			dataframe = dataframe.rename(columns={'gene_id': 'Gene ID', 'gene_name': 'Gene name', 'baseMean': 'Average expression', 'pvalue': 'P-value', 'padj': 'Adjusted P-value', 'differential_expression': 'Differentially expressed'}).drop(['significant'], axis=1)

			# Write
			dataframe.to_excel(writer, sheet_name=comparison, index=False)
			
			# Modify
			workbook = writer.book
			worksheet = writer.sheets[comparison]
			format1 = workbook.add_format({'align': 'center'})
			worksheet.set_column('A:J', 20, format1)
			max_logfc = dataframe['log2FoldChange'].abs().max()
			worksheet.conditional_format('D1:D'+str(len(dataframe.index)+1), {'type': '3_color_scale', 'min_value': -max_logfc, 'min_type': 'num', 'mid_value': 0, 'mid_type': 'num', 'max_value': max_logfc, 'max_type': 'num', 'min_color': '#67a9cf', 'mid_color': '#ffffff', 'max_color': '#ef8a62'})
			worksheet.autofilter(0, 0, len(dataframe.index), len(dataframe.columns) - 1)

#######################################################
#######################################################
########## S6. Enrichment
#######################################################
#######################################################

#############################################
########## 1. GO
#############################################

# @follows(runDESeq2)

@transform('arion/illumina/s05-differential_expression.dir/*/*-gene-deseq.tsv',
		   regex(r'(.*)/s05-differential_expression.dir/(.*)/(.*)-gene-deseq.tsv'),
		   add_inputs(r'arion/isoseq/s05-talon.dir/\2/gtf/*-SJ_filtered.gtf'),
		   r'\1/s06-enrichment.dir/\2/go/\3-go_enrichment.tsv')

def runGoEnrichment(infiles, outfile):

	# Run
	run_r_job('run_go_enrichment', infiles, outfile, conda_env='env', W='00:30', GB=50, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 2. Novelty
#############################################

# @follows(runDESeq2)

# @transform('arion/illumina/s05-differential_expression.dir/*/*-gene-deseq.tsv',
# 		   regex(r'(.*)/s05-differential_expression.dir/(.*)/(.*)-gene-deseq.tsv'),
# 		   add_inputs(r'arion/isoseq/s05-talon.dir/\2/*abundance_filtered.tsv'),
# 		   r'\1/s06-enrichment.dir/\2/novelty/\3-novelty_enrichment.rda')

# def runNoveltyEnrichment(infiles, outfile):

# 	# Run
# 	run_r_job('run_novelty_enrichment', infiles, outfile, modules=['R/4.0.3'], W='01:00', GB=15, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

# rm -r arion/illumina/s06-enrichment.dir/human/novelty

#############################################
########## 3. Transcript Novelty
#############################################

# @follows(runDESeq2)

@transform('arion/illumina/s05-differential_expression.dir/*human*/*-transcript-deseq.tsv',
		   regex(r'(.*)/s05-differential_expression.dir/(.*)/(.*)-transcript-deseq.tsv'),
		   add_inputs(r'arion/isoseq/s05-talon.dir/\2/gtf/*-SJ_filtered-transcript_classification.tsv'),
		   r'\1/s06-enrichment.dir/\2/transcript_novelty/\3-transcript_novelty_enrichment.rda')

def runTranscriptNoveltyEnrichment(infiles, outfile):

	# Run
	run_r_job('run_transcript_novelty_enrichment', infiles, outfile, modules=['R/4.0.3'], W='30:00', GB=30, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 4. Domain
#############################################

# @follows(runDESeq2)

@transform('arion/illumina/s05-differential_expression.dir/*/*-transcript-deseq.tsv',
		   regex(r'(.*)/s05-differential_expression.dir/(.*)/(.*)-transcript-deseq.tsv'),
		   add_inputs(r'arion/isoseq/s07-pfam.dir/\2/\2-translated_pfam.tsv'),
		   r'\1/s06-enrichment.dir/\2/domain/\3-domain_enrichment.tsv')

def runDomainEnrichment(infiles, outfile):

	# Run
	run_r_job('run_domain_enrichment', infiles, outfile, conda_env='env', W='00:15', GB=15, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

# #############################################
# ########## 5. Repeats
# #############################################

# @follows(runDESeq2)

@transform('arion/illumina/s05-differential_expression.dir/*/*-transcript-deseq.tsv',
		   regex(r'(.*)/s05-differential_expression.dir/(.*)/(.*)-transcript-deseq.tsv'),
		   add_inputs(r'arion/isoseq/s08-repeatmasker.dir/\2/*_repeatmasker.tsv'),
		   r'\1/s06-enrichment.dir/\2/repeat/\3-repeat_enrichment.tsv')

def runRepeatEnrichment(infiles, outfile):

	# Run
	run_r_job('run_repeat_enrichment', infiles, outfile, conda_env='env', W='00:30', GB=15, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 5. miRNA
#############################################

# @follows(runDESeq2)

@transform('arion/illumina/s05-differential_expression.dir/human/*-transcript-deseq.tsv',
		   regex(r'(.*)/s05-differential_expression.dir/(.*)/(.*)-transcript-deseq.tsv'),
		   add_inputs(r'arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance_scan/mirna/miranda-hits_filtered_80pct.tsv'),
		   r'\1/s06-enrichment.dir/\2/mirna/\3-mirna_enrichment_80pct.tsv')

def runMirnaEnrichment(infiles, outfile):

	# Run
	run_r_job('run_mirna_enrichment', infiles, outfile, conda_env='env', W='00:30', GB=15, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 5. miRNA
#############################################

# @follows(runDESeq2)

@transform('arion/illumina/s05-differential_expression.dir/human/*-gene-deseq.tsv',
		   regex(r'(.*)/s05-differential_expression.dir/(.*)/(.*)-gene-deseq.tsv'),
		   add_inputs(r'arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance_scan/mirna/miranda-hits_filtered_85pct.tsv'),
		   r'\1/s06-enrichment.dir/\2/mirna_geneid_default/\3-mirna_enrichment_85pct.tsv')

def runMirnaEnrichmentGene(infiles, outfile):

	# Run
	run_r_job('run_mirna_enrichment_gene', infiles, outfile, conda_env='env', W='01:00', GB=50, n=1, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#######################################################
#######################################################
########## S7. WGCNA
#######################################################
#######################################################

#############################################
########## 1. Pick soft thresholds
#############################################

# @follows(getGeneExpression)

@transform('arion/illumina/s04-expression.dir/human/human-counts.rda',
		   regex(r'(.*)/s04-expression.dir/(.*)/.*.rda'),
		   r'\1/s07-wgcna.dir/\2/network/\2-soft_thresholds_signed.rda')

def pickSoftThresholds(infile, outfile):

	# Run
	run_r_job('pick_soft_thresholds', infile, outfile, modules=['R/4.0.3'], W='00:45', GB=25, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 2. Cluster genes
#############################################

# @follows(getGeneExpression)

@transform(pickSoftThresholds,
		   regex(r'(.*)-soft_thresholds_(.*).rda'),
		   r'\1-gene_network_\2.rda')

def clusterGenes(infile, outfile):

	# Run
	run_r_job('cluster_genes', infile, outfile, modules=['R/4.0.3'], W='00:45', GB=50, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 3. Get modules
#############################################

# @follows(getGeneExpression)

@transform(clusterGenes,
		   suffix('.rda'),
		   add_inputs(pickSoftThresholds),
		   '_modules.rda')

def getGeneModules(infiles, outfile):

	# Run
	run_r_job('get_gene_modules', infiles, outfile, modules=['R/4.0.3'], W='00:45', GB=15, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 4. Get enrichment
#############################################

# @follows(getGeneExpression)

@transform(getGeneModules,
		   suffix('s.rda'),
		   '_enrichment.rda')

def runModuleEnrichment(infiles, outfile):

	# Run
	run_r_job('run_module_enrichment', infiles, outfile, run_locally=True)# conda env, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 5. Get module preservation
#############################################

# @follows(getGeneExpression)

@transform(('arion/geo_illumina/s04-expression.dir/*/*-counts.rda', 'arion/geo_illumina/s05-primates.dir/*/RSEM/counts/*-counts.rda'),
		   regex(r'.*/(.*)-counts.rda'),
		   add_inputs(pickSoftThresholds, getGeneModules),
		   r'arion/illumina/s07-wgcna.dir/human/module_preservation/\1-module_preservation.rda')

def getModulePreservation(infiles, outfile):

	# Run
	run_r_job('get_module_preservation', infiles, outfile, W='02:00', GB=10, n=5, modules=['R/4.0.3'], stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 5. Get novel gene signatures
#############################################

# @follows(getGeneExpression)

@transform('arion/illumina/s07-wgcna.dir/human/network/human-gene_network_signed_adjacency.rda',
		   regex(r'.*/(.*)-gene_network_signed_adjacency.rda'),
		   add_inputs(clusterNovelGenes),
		   r'arion/illumina/s07-wgcna.dir/human/network/\1-novel_gene_signatures.rda')

def getNovelGeneSignatures(infiles, outfile):

	# Run
	run_r_job('get_novel_gene_signatures', infiles, outfile, W='02:00', GB=50, n=1, modules=['R/4.0.3'], stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 6. Get gene connectivity
#############################################

# @follows(getGeneExpression)

# @transform(getGeneModules,
# 		   suffix('s.rda'),
# 		   add_inputs(pickSoftThresholds, 'arion/illumina/s07-wgcna.dir/human/network/human-gene_network_signed_adjacency.rda'),
# 		   '_membership.rda')

# def getModuleMembership(infiles, outfile):

# 	# Run
# 	run_r_job('get_module_membership', infiles, outfile, W='00:10', GB=30, n=1, modules=['R/4.0.3'], stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 6. Filter adjacency
#############################################

@subdivide('arion/illumina/s07-wgcna.dir/human/network/human-gene_network_signed_adjacency.rda',
		   regex(r'(.*)/network/(.*).rda'),
		   r'\1/graph/networks/\2_*.rda',
		   r'\1/graph/networks/\2_{min_adjacency}.rda')

def filterAdjacency(infile, outfiles, outfileRoot):

	# Loop
	for min_adjacency in [.07, 0.071, 0.072, 0.073, 0.074, 0.075, 0.076, 0.077, 0.078, 0.079, 0.08, 0.081, 0.082, 0.083, 0.084, 0.085, 0.086, 0.087, 0.088, 0.089, 0.09, 0.091, 0.092, 0.093, 0.094, 0.095, 0.096, 0.097, 0.098, 0.099, 0.1]:

		# Get outfile
		outfile = outfileRoot.format(**locals())

		# Run
		run_r_job('filter_adjacency', infile, outfile, q='sla', additional_params=min_adjacency, modules=['R/4.0.3'], W='10:00', GB=10, n=10, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

# find arion/illumina/s07-wgcna.dir/human/graph/networks -name "*adjacency*.log" | js

#############################################
########## 7. Filter adjacency
#############################################

@subdivide('arion/illumina/s07-wgcna.dir/human/graph/networks/*08*.rda',
		   regex(r'(.*)/.*/(.*).rda'),
		   add_inputs('arion/illumina/s07-wgcna.dir/human/network/human-gene_network_signed_modules.rda'),
		   r'\1/plots/\2_*.tsv.gz',
		   r'\1/plots/\2_{layout}_{network_type}_{random_seed}.tsv.gz')

def plotNetwork(infiles, outfiles, outfileRoot):

	# Loop
	for random_seed in [x for x in range(10)]:
		for network_type in ['full']:
			for layout in ['fruchtermanreingold']:

				# Get outfile
				outfile = outfileRoot.format(**locals())

				# Run
				run_r_job('plot_network', infiles, outfile, run_locally=False, q='sla', print_outfile=False, additional_params=[random_seed, network_type, layout], modules=['R/4.0.3'], W='03:00', GB=20, n=5, stdout=outfile.replace('.tsv.gz', '.log'), stderr=outfile.replace('.tsv.gz', '.err'))

#############################################
########## 8. Plot
#############################################

@subdivide('arion/illumina/s07-wgcna.dir/human/graph/old/plots-4/human-gene_network_signed_adjacency_0.1_fruchtermanreingold_full_64.tsv.gz',
		   regex(r'(.*)/old.*/(.*).tsv.gz'),
		   add_inputs('arion/illumina/s07-wgcna.dir/human/network/human-gene_network_signed_modules.rda'),
		   r'\1/plots_final/\2_*.tsv.gz',
		   r'\1/plots_final/\2_size_{edge_size}_edgealpha_{edge_alpha}_curvature_{curvature}.png')

def plotGraph(infiles, outfiles, outfileRoot):

	# Loop
	for edge_size in [0.5]:
		for edge_alpha in [0.005]:
			for curvature in [-0.1, 0, 0.1]:
			# for curvature in [round(x, 2) for x in np.arange(-0.1, 0.1, 0.01)]:

				# Get outfile
				outfile = outfileRoot.format(**locals())

				# Run
				run_r_job('plot_graph', infiles, outfile, run_locally=False, q='sla', print_outfile=False, additional_params=[edge_size, edge_alpha, curvature], modules=['R/4.0.3'], W='01:00', GB=10, n=5)#, stdout=outfile.replace('.tsv.gz', '.log'), stderr=outfile.replace('.tsv.gz', '.err'))

#############################################
########## 8. Get module TSSs
#############################################

@subdivide(getGeneModules,
		   regex(r'(.*)/network/(.*)-.*rda'),
		   add_inputs(r'arion/isoseq/s05-talon.dir/\2/gtf/*_talon-SJ_filtered.gtf', r'arion/atacseq/s05-counts.dir/\2/\2-atacseq_consensus_peaks_3kb.bed'),
		   r'\1/motifs/peaks/tss_intersect_3kb/\2-module_*.bed',
		   r'\1/motifs/peaks/tss_intersect_3kb/\2-{selected_module}-peaks.bed')

def getModulePeaks(infiles, outfiles, outfileRoot):

	# Run
	run_r_job('get_module_peaks', infiles, outfileRoot, run_locally=False, conda_env='env', W='00:10', GB=10, n=1, modules=['R/4.0.3'], stdout=os.path.join(os.path.dirname(outfileRoot), 'job.log'), stderr=os.path.join(os.path.dirname(outfileRoot), 'job.err'))

#############################################
########## 9. Get module motifs
#############################################

# @transform(getModulePeaks,
@transform('arion/illumina/s07-wgcna.dir/human/motifs/peaks/*/*module*.bed',
		   regex(r'(.*)/peaks/(.*?)(-.*)-peaks.bed'),
		   add_inputs(r'arion/datasets/reference_genomes/human/*.dna_sm.primary_assembly.fa', r'\1/peaks/\2*background*.bed'),
		   r'\1/homer_given/\2\3/knownResults.html')

def getModuleMotifs(infiles, outfile):

	# Command (add genome)
	outdir = os.path.dirname(outfile)
	cmd_str = ''' findMotifsGenome.pl {infiles[0]} {infiles[1]} {outdir}/ -size given -bg {infiles[2]} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, print_cmd=False, modules=['homer/4.10'], W='06:00', GB=10, n=1, jobname=outfile.split('/')[-2]+'-homer', stdout=outfile.replace('knownResults.html', 'job.log'), stderr=outfile.replace('knownResults.html', 'job.err'))

#######################################################
#######################################################
########## S9. SUPPA
#######################################################
#######################################################

#############################################
########## 1. Index
#############################################

def suppaIndexJobs():
	for organism, organism_references in reference_dict.items():
		for file_format in ['ioi', 'ioe']:
			infile = organism_references['isoseq']['gtf_filtered']
			outdir = 'arion/illumina/s08-suppa.dir/{organism}/01-events/{file_format}/'.format(**locals())
			yield [infile, outdir, file_format]

# @follows(functionToFollow)

@files(suppaIndexJobs)

def buildSuppaIndex(infile, outdir, file_format):

	# Basename
	basename = outdir.split('/')[-4]

	# Command
	if file_format == 'ioe':
		cmd_str = '''python $SUPPA_HOME/suppa.py generateEvents -i {infile} -o {outdir}{basename} -f {file_format} -e SE SS MX RI FL'''.format(**locals())
	elif file_format == 'ioi':
		cmd_str = '''python $SUPPA_HOME/suppa.py generateEvents -i {infile} -o {outdir}{basename} -f {file_format}'''.format(**locals())

	# Run
	run_job(cmd_str, outdir, W='00:15', modules=['suppa/2.3'], GB=10, n=1, stdout=os.path.join(outdir, 'job.log'), jobname='_'.join(outdir.split('/')[-4:]).replace('_01-events', ''), ow=True)

#############################################
########## 2. PSI
#############################################

def psiJobs():
	for organism in ['human', 'mouse']:
		tpm_file = 'arion/illumina/s04-expression.dir/{organism}/{organism}-transcript_tpm.txt'.format(**locals())
		indices = glob.glob('arion/illumina/s08-suppa.dir/{organism}/01-events/io*/*.io?'.format(**locals()))
		for suppa_index in indices:
			file_format = suppa_index.split('.')[-1]
			if file_format == 'ioe':
				infiles = [tpm_file, suppa_index]
				event_type = suppa_index.split('_')[-2]
			elif file_format == 'ioi':
				infiles = [tpm_file, reference_dict[organism]['isoseq']['gtf_filtered']]
				event_type = 'isoform'
			outfile = 'arion/illumina/s08-suppa.dir/{organism}/02-psi/{organism}_{event_type}.psi'.format(**locals())
			yield [infiles, outfile, event_type] 

@follows(buildSuppaIndex)

@files(psiJobs)

def getSuppaPSI(infiles, outfile, event_type):
	
	# Isoform
	if event_type == 'isoform':
		outname = outfile[:-len('_isoform.psi')]
		cmd_str = 'python $SUPPA_HOME/suppa.py psiPerIsoform -g {infiles[1]} -e {infiles[0]} -o {outname}'.format(**locals())
	# Event
	else:
		outname = outfile[:-len('.psi')]
		cmd_str = 'python $SUPPA_HOME/suppa.py psiPerEvent --ioe-file {infiles[1]} --expression-file {infiles[0]} -o {outname}'.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="00:30", GB=10, n=1, modules=['suppa/2.3'], ow=False, stderr=outfile.replace('.psi', '.err'))

#############################################
########## 3. Split
#############################################

def splitJobs():
	infiles = glob.glob('arion/illumina/s04-expression.dir/*/*-transcript_tpm.txt')+glob.glob('arion/illumina/s08-suppa.dir/*/02-psi/*.psi')
	for infile in infiles:
		infile_split = infile.split('/')
		organism = infile_split[3]
		basename = infile_split[-1].replace('-', '_').replace('.psi', '').replace('.txt', '')#.replace('all_', '')
		outdir = 'arion/illumina/s08-suppa.dir/{organism}/03-split/{basename}'.format(**locals())
		yield [infile, outdir, basename]

# @follows(getTranscriptTPM, getSuppaPSI)

@files(splitJobs)

def splitData(infile, outdir, basename):

	# Get metadata
	dataframe = pd.read_table(infile)

	# Get groups
	group_dataframe = pd.DataFrame(dataframe.columns).rename(columns={0: 'sample'})
	group_dataframe['cell_type'] = [x.split('_')[1].replace('2PN', '1C') for x in group_dataframe['sample']]
	group_dict = group_dataframe.groupby('cell_type')['sample'].apply(lambda x: list(x))

	# Loop
	for cell_type, samples in group_dict.items():

		# Get outfile
		file_format = 'tsv' if 'tpm' in infile else 'psi'
		outfile = '{outdir}/{basename}_{cell_type}.{file_format}'.format(**locals())

		# Outdir
		if not os.path.exists(outdir):
			os.makedirs(outdir)
		
		# Subset and write
		dataframe[samples].to_csv(outfile, sep='\t', index_label=False, na_rep='nan')

#############################################
########## 4. Differential splicing
#############################################

def suppaRunJobs():
	for organism, comparisons in comparison_dict.items():
		for event_type in ['isoform', 'A3', 'A5', 'AF', 'AL', 'MX', 'RI', 'SE']:
			if event_type == 'isoform':
				suppa_index = 'arion/illumina/s08-suppa.dir/{organism}/01-events/ioi/{organism}.ioi'.format(**locals())
			else:
				suppa_index = 'arion/illumina/s08-suppa.dir/{organism}/01-events/ioe/{organism}_{event_type}_strict.ioe'.format(**locals())
			for comparison in comparisons:
				infiles = []
				for cell_type in comparison:
					infiles.append('arion/illumina/s08-suppa.dir/{organism}/03-split/{organism}_{event_type}/{organism}_{event_type}_{cell_type}.psi'.format(**locals()))
					infiles.append('arion/illumina/s08-suppa.dir/{organism}/03-split/{organism}_transcript_tpm/{organism}_transcript_tpm_{cell_type}.tsv'.format(**locals()))
				infiles.append(suppa_index)
				outfile = 'arion/illumina/s08-suppa.dir/{organism}/04-dpsi/{organism}-{comparison[0]}_vs_{comparison[1]}/{organism}-{comparison[0]}_vs_{comparison[1]}-{event_type}/{organism}-{comparison[0]}_vs_{comparison[1]}-{event_type}.psivec'.format(**locals())
				yield [infiles, outfile]

# @follows(splitData)

@files(suppaRunJobs)

def getDifferentialPSI(infiles, outfile):

	# Full paths
	infiles = [os.path.join(os.getcwd(), x) for x in infiles]

	# Command
	outname = os.path.join(os.getcwd(), outfile.rsplit('.', 1)[0])
	
	# Command
	cmd_str = 'python $SUPPA_HOME/suppa.py diffSplice -m empirical --save_tpm_events --input {infiles[4]} --psi {infiles[0]} {infiles[2]} --tpm {infiles[1]} {infiles[3]} -gc -o {outname}'.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="06:00", GB=50, n=1, modules=['suppa/2.3', 'python/3.8.2'], stdout=outfile.replace('.psivec', '.log'), stderr=outfile.replace('.psivec', '.err'))

# find arion/illumina/s08-suppa.dir/*/04-dpsi -name "*.log" | jsc
# find arion/illumina/s08-suppa.dir/*/04-dpsi -name "*.log" | xargs wc -l
# find arion/illumina/s08-suppa.dir/*/04-dpsi -name "*.err" | xargs wc -l

# # ls /hpc/users/torred23/pipelines/projects/early-embryo/arion/illumina/s08-suppa.dir/*/*/04-dpsi/*/*/*.log | js | grep -v completed
# find arion/illumina/s08-suppa.dir/mouse/ensembl/04-dpsi -name "*.psivec" | wc -l
# # ls /hpc/users/torred23/pipelines/projects/early-embryo/arion/illumina/s08-suppa.dir/*/*/04-dpsi/*/*/*.err | lr | grep Error > error_samples.txt

#############################################
########## 5. Summary
#############################################

# getDifferentialPSI
@collate('arion/illumina/s08-suppa.dir/*/04-dpsi/*/*/*.psivec',
		 regex(r'(.*)/04-dpsi/(.*)/.*/.*.psivec'),
		 r'\1/05-summaries/\2-suppa_summary.tsv')

def createSuppaSummary(infiles, outfile):

	# Make summary
	summary_dataframe = pd.concat([S.createSuppaSummary(x[:-len('.psivec')]) for x in infiles if 'isoform' not in x])

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Write
	summary_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 6. Summary
#############################################

# getDifferentialPSI
@transform('arion/illumina/s08-suppa.dir/*/04-dpsi/*/*/*isoform.psivec',
		    regex(r'(.*)/04-dpsi/(.*)/.*/.*.psivec'),
		    r'\1/05-summaries/\2-suppa_summary_isoform.tsv')

def createSuppaIsoformSummary(infile, outfile):

	# Make summary
	summary_dataframe = S.createSuppaSummary(infile[:-len('.psivec')], event_type='isoform')

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Write
	summary_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 6. Cluster PSI
#############################################

@collate(createSuppaIsoformSummary,
		 regex(r'(.*)/05-summaries/(human)-.*.tsv'),
		 add_inputs(r'\1/02-psi/\2_isoform.psi'),
		 r'\1/06-psi_clusters_isoswitch_significance/\2_isoform-timepoints.rda')

def clusterPSI(infiles, outfile):

	# Split
	summary_files = [x[0] for x in infiles]
	psi_file = infiles[0][1]
	infiles = [psi_file]+summary_files

	# Run
	run_r_job('cluster_psi', infiles, outfile, modules=['R/4.0.3'], W="02:00", GB=15, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 7. Get clusters
#############################################

@transform(clusterPSI,
		   suffix('timepoints.rda'),
		   'clusters.rda')

def getPSIClusters(infile, outfile):

	# Run
	run_r_job('get_psi_clusters', infile, outfile, modules=['R/4.0.3'], W="01:00", GB=15, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 8. Correlate gene expression
#############################################

@transform('arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance/human_isoform-clusters.rda',
		   regex(r'(.*)/(.*)-clusters.rda'),
		   add_inputs('arion/illumina/s08-suppa.dir/human/02-psi/human_isoform.psi', 'arion/illumina/s04-expression.dir/human/human-gene_normalized_counts.tsv'),
		   r'\1/gene_correlation/\2_cluster_gene_correlation.tsv')

def getPSIClusterCorrelations(infiles, outfile):

	# Run
	run_r_job('get_psi_cluster_correlations', infiles, outfile, modules=['R/4.0.3'], W='10:00', GB=30, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 8. Correlate gene expression to SE
#############################################

@transform('arion/illumina/s08-suppa.dir/human/02-psi/human_SE.psi',
		   regex(r'(.*)/02-psi/(.*).psi'),
		   add_inputs('arion/illumina/s04-expression.dir/human/human-gene_normalized_counts.tsv', 'arion/illumina/s08-suppa.dir/human/05-summaries/human-*-suppa_summary.tsv'),
		   r'\1/07-se_psi_gene_correlation/\2_gene_correlation_rbp.tsv')

def getAsExonGeneCorrelations(infiles, outfile):

	# Run
	run_r_job('get_as_exon_gene_correlations', infiles, outfile, modules=['R/4.0.3'], W='01:00', GB=30, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 9. Get cluster info
#############################################

@subdivide('arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance/human_isoform-clusters.rda',
		   regex(r'(.*)/(.*)/(.*)-clusters.rda'),
		   add_inputs('arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.gtf', 'arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.fasta'),
		   r'\1/\2_scan/*/*',
		   r'\1/\2_scan/cluster_data/\3-{cluster}.{format}')

def getPSIClusterData(infiles, outfiles, outfileRoot):

	# Run
	run_r_job('get_psi_cluster_data', infiles, outfileRoot, modules=['R/4.0.3'], W='01:00', GB=30, n=1)

#############################################
########## 10. Scan RBP overlap
#############################################

@transform('arion/datasets/encode/eclip/links/*.bed.gz',
		   regex(r'.*/(.*).bed.gz'),
		   add_inputs('arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance_scan/cluster_data/human_isoform-cluster_*.bed'),
		   r'arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance_scan/rbp/\1-isoform_cluster_overlap.bed')

def scanRbpOverlap(infiles, outfile):

	# Run
	run_r_job('scan_rbp_overlap', infiles, outfile, modules=['R/4.0.3'], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

# find /hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/encode/eclip/links -name "*POL2RG*"

#############################################
########## 10. Download
#############################################

@originate('arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance_scan/mirna/fasta/mature.fa')

def downloadMirnaFasta(outfile):

	# Get outdir
	outdir = os.path.dirname(outfile)

	# Command
	cmd_str = ''' cd {outdir} && wget --no-check-certificate https://www.mirbase.org/ftp/CURRENT/mature.fa.gz && gunzip mature.fa.gz '''.format(**locals())

	# Run
	# run_job(cmd_str, outfile, run_locally=True)

#############################################
########## 10. Subset
#############################################

@transform(downloadMirnaFasta,
		   suffix('.fa'),
		   '.hsa.fa')

def subsetHumanMirna(infile, outfile):

	# Run
	pass
	# run_r_job('subset_human_mirna', infile, outfile, modules=['R/4.0.3'], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 10. Run miRanda
#############################################

# @transform('arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance_scan/cluster_data/human_isoform-*.fasta',
# 		   regex(r'(.*)/cluster_data/(.*)-(.*).fasta'),
# 		   add_inputs(subsetHumanMirna),
# 		   r'\1/mirna/results/\3/\2-\3_miranda.txt')

@transform('arion/isoseq/s08-repeatmasker.dir/human/fasta/Homo_sapiens.GRCh38.102_talon-SJ_filtered.fasta.*',
		   regex(r'.*/(.*)'),
		   add_inputs(subsetHumanMirna),
		   r'arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance_scan/mirna/results_all/\1/\1_miranda.txt')

def runMiranda(infiles, outfile):

	# Parameters
	hit_file = outfile.replace('.txt', '.hit.txt')
	scan_file = outfile.replace('.txt', '.scan.txt')

	# Command
	cmd_str = ''' miranda {infiles[1]} {infiles[0]} -out {outfile} && \
		grep -e '^>hsa' {outfile} > {hit_file} && \
		grep -e '^>>hsa' {outfile} > {scan_file}
	'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['miranda/3.3a'], W='06:00', GB=15, n=1, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'), print_outfile=False)

#############################################
########## 10. Filter miRanda
#############################################

@transform('arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance_scan/mirna/results_all/*/*.hit.txt',
		   regex(r'(.*).txt'),
		   add_inputs('arion/isoseq/s06-cpat.dir/human/Homo_sapiens.GRCh38.102_talon-SJ_filtered-cpat.ORF_prob.best.tsv', 'arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.gtf'),
		   r'\1_filtered_geneid_80pct.tsv')

def filterMirandaResults(infiles, outfile):

	# Run
	run_r_job('filter_miranda_results', infiles, outfile, conda_env='env', W='00:15', GB=15, n=1, run_locally=False)#, print_outfile=False, print_cmd=False)

#############################################
########## 10. Create geneset
#############################################

@collate('arion/illumina/s08-suppa.dir/human/06-psi_clusters_suppa_significance_scan/mirna/results_all/*/*.hit_filtered_geneid.tsv',
		 regex(r'(.*)/results_all/.*.tsv'),
		 r'\1/miranda-hits_filtered_85pct.tsv')

def filterMirandaHits(infiles, outfile):

	# Read
	result_dataframe = pd.concat([pd.read_table(x).query('pct_identity_wobble>85') for x in infiles])[['mirna', 'gene_id', 'transcript_id']].drop_duplicates()

	# Writw
	result_dataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S10. Isoform switching
#######################################################
#######################################################

#############################################
########## 1. Filter
#############################################

def isoformFilterJobs():
	for organism in ['human']: #, 'mouse'
		filtered_gtf = glob.glob('arion/isoseq/s05-talon.dir/{organism}/gtf/*.102_talon-SJ_filtered.gtf'.format(**locals()))[0]
		for file_type in ['gtf_cds', 'cpat_predictions', 'pfam_predictions', 'transcript_fasta']:
			infile = reference_dict[organism]['isoseq'][file_type]
			infile_basename, infile_extension = os.path.basename(infile).rsplit('.', 1)
			infiles = [infile, filtered_gtf]
			outfile = 'arion/illumina/s09-isoform_switching.dir/{organism}/filtered_data/{infile_basename}_filtered.{infile_extension}'.format(**locals())
			yield [infiles, outfile, file_type]

@files(isoformFilterJobs)

def filterIsoformData(infiles, outfile, file_type):

	# Run
	run_r_job('filter_isoform_data', infiles, outfile, additional_params=file_type, conda_env='env')

#############################################
########## 2. Load data
#############################################

@collate('arion/illumina/s09-isoform_switching.dir/human/filtered_data/*',
		 regex(r'(.*)/(.*)/filtered_data/.*'),
		 add_inputs(r'arion/illumina/s04-expression.dir/\2/\2-sample_metadata.txt', comparison_file),
		 r'\1/\2/\2-isoforms.rda')

def loadIsoformData(infiles, outfile):

	# Split
	infiles = [x[0] for x in infiles]+[infiles[0][1], infiles[0][2]]

	# Run
	run_r_job('load_isoform_data', infiles, outfile, conda_env='env', W="06:00", GB=30, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 3. Run
#############################################

@transform(loadIsoformData,
		  suffix('.rda'),
		  '_results.rda')

def getIsoformSwitching(infile, outfile):

	# Run
	run_r_job('get_isoform_switching', infile, outfile, conda_env='env', W="06:00", GB=30, n=1, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 4. Cluster PSI IsoSwitch
#############################################

@transform(getIsoformSwitching,
		   suffix('.rda'),
		   '_psi_clusters.rda')

def clusterPsiIsoSwitch(infile, outfile):

	# Run
	run_r_job('cluster_psi_isoswitch', infile, outfile, modules=['R/4.0.3'], W="02:00", GB=15, n=1, run_locally=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#######################################################
#######################################################
########## S10. Promoter activity
#######################################################
#######################################################

#############################################
########## 1. Run
#############################################

@collate('arion/illumina/s03-alignment.dir/human/isoseq/STAR/pass2/*/*-SJ.out.tab',
		 regex(r'(.*)/s03-alignment.dir/(.*?)/.*.tab'),
		 add_inputs(r'arion/isoseq/s05-talon.dir/\2/gtf/*102_talon-SJ_filtered.gtf'),
		 r'\1/s10-promoter_activity.dir/\2/\2-promoter_activity_grouped.rda')

def runProActiv(infiles, outfile):

	# Split infiles
	gtf_file = infiles[0][1]
	infiles = [x[0] for x in infiles]

	# Run
	run_r_job('run_proactiv', infiles, outfile, additional_params=gtf_file, modules=['R/4.0.3'], W='01:00', GB=10, n=5, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))


#######################################################
#######################################################
########## S11. RBP scan
#######################################################
#######################################################

#############################################
########## 1. Download motifs
#############################################

@originate('arion/illumina/s11-motif_scan.dir/motifs/ATtRACT/pwm.txt')

def downloadAttractMotifs(outfile):

	# Command
	outdir = os.path.dirname(outfile)
	cmd_str = ''' mkdir -p {outdir} && cd {outdir} && wget https://attract.cnic.es/attract/static/ATtRACT.zip && unzip ATtRACT.zip '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, run_locally=True)

#############################################
########## 2. Split
#############################################

@subdivide(downloadAttractMotifs,
		   regex(r'(.*)/pwm.txt'),
		   add_inputs(r'\1/ATtRACT_db.txt'),
		   r'\1/*/*/*.motif',
		   r'\1/split/{motif_name}_{organism}.motif')

def splitAttractMotifs(infiles, outfiles, outfileRoot):

	# Run
	run_r_job('split_attract_motifs', infiles, outfileRoot, run_locally=True, print_outfile=False)#, modules=['R/4.0.3'], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 3. Scan motifs
#############################################

# @transform(splitAttractMotifs,
@transform('/hpc/users/torred23/pipelines/projects/early-embryo/embryo-illumina/arion/illumina/s11-motif_scan.dir/motifs/ATtRACT/split/*human.motif',
		   regex(r'(.*)/motifs/.*/(.*)_(.*).motif'),
		   add_inputs(r'/hpc/users/torred23/pipelines/projects/early-embryo/embryo-illumina/arion/datasets/reference_genomes/\3/*.dna_sm.primary_assembly.fa'),
		   r'\1/\3/homer/\2_\3/\2_\3_homer_results.txt')

def findRbpMotifs(infiles, outfile):

	# Command
	outdir = os.path.dirname(outfile)
	cmd_str = ''' cd {outdir} && scanMotifGenomeWide.pl {infiles[0]} {infiles[1]} -keepAll > {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['homer/4.10'], W='03:00', GB=50, n=1, print_outfile=True, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'))

# find arion/illumina/s11-motif_scan.dir/human/homer -name "*.log" | js | grep -v completed > failed_homer_scan.sh
# find arion/illumina/s11-motif_scan.dir/human/homer -name "*.log" | wc -l

#############################################
########## 4. Get isoform BED
#############################################

@transform('arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.gtf',
		   regex(r'.*/(.*)/gtf/(.*).gtf'),
		   r'arion/illumina/s11-motif_scan.dir/\1/bed/\2_isoform.bed')

def getIsoformBed(infile, outfile):

	# Run
	run_r_job('get_isoform_bed', infile, outfile, modules=['R/4.0.3'], W='00:15', GB=10, n=1, run_locally=True, print_outfile=False, print_cmd=False)

#############################################
########## 5. Intersect isoforms and motifs
#############################################

@transform(('arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000015479*/*ENSG00000015479*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000092199*/*ENSG00000092199*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000100650*/*ENSG00000100650*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000100883*/*ENSG00000100883*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000104824*/*ENSG00000104824*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000112081*/*ENSG00000112081*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000116001*/*ENSG00000116001*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000128016*/*ENSG00000128016*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000136527*/*ENSG00000136527*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000140319*/*ENSG00000140319*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000141543*/*ENSG00000141543*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000147274*/*ENSG00000147274*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000152795*/*ENSG00000152795*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000153187*/*ENSG00000153187*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000161547*/*ENSG00000161547*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000188529*/*ENSG00000188529*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000265241*/*ENSG00000265241*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000020577*/*ENSG00000020577*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000048740*/*ENSG00000048740*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000065978*/*ENSG00000065978*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000070756*/*ENSG00000070756*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000078328*/*ENSG00000078328*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000080802*/*ENSG00000080802*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000084072*/*ENSG00000084072*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000089280*/*ENSG00000089280*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000096746*/*ENSG00000096746*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000099783*/*ENSG00000099783*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000102317*/*ENSG00000102317*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000107105*/*ENSG00000107105*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000107201*/*ENSG00000107201*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000111786*/*ENSG00000111786*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000114416*/*ENSG00000114416*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000115875*/*ENSG00000115875*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000116560*/*ENSG00000116560*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000119707*/*ENSG00000119707*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000121774*/*ENSG00000121774*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000122566*/*ENSG00000122566*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000124193*/*ENSG00000124193*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000125207*/*ENSG00000125207*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000131503*/*ENSG00000131503*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000131914*/*ENSG00000131914*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000132485*/*ENSG00000132485*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000135316*/*ENSG00000135316*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000135829*/*ENSG00000135829*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000136231*/*ENSG00000136231*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000136450*/*ENSG00000136450*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000138668*/*ENSG00000138668*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000138757*/*ENSG00000138757*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000143889*/*ENSG00000143889*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000144642*/*ENSG00000144642*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000147140*/*ENSG00000147140*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000148584*/*ENSG00000148584*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000151923*/*ENSG00000151923*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000151962*/*ENSG00000151962*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000153037*/*ENSG00000153037*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000160710*/*ENSG00000160710*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000164548*/*ENSG00000164548*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000164902*/*ENSG00000164902*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000165119*/*ENSG00000165119*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000168066*/*ENSG00000168066*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000169045*/*ENSG00000169045*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000170144*/*ENSG00000170144*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000197111*/*ENSG00000197111*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000197451*/*ENSG00000197451*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000204356*/*ENSG00000204356*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000214575*/*ENSG00000214575*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000239306*/*ENSG00000239306*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000135486*/*ENSG00000135486*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000159217*/*ENSG00000159217*_homer_results.txt', 'arion/illumina/s11-motif_scan.dir/human/homer/*ENSG00000169813*/*ENSG00000169813*_homer_results.txt'),
		   regex(r'(.*)/homer/(.*)/.*.txt'),
		   add_inputs(getIsoformBed),
		   r'\1/intersect/\2/\2-isoform_intersect_coordinates.tsv')

def intersectIsoformMotifs(infiles, outfile):

	# Run
	run_r_job('intersect_isoform_motifs', infiles, outfile, modules=['R/4.0.3'], W='00:15', GB=15, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 6. Collapse RBPs
#############################################

@collate('arion/illumina/s11-motif_scan.dir/human/intersect/*/*-isoform_intersect.tsv',
		 regex(r'(.*)/intersect/(.*?_.*?)_.*/.*-(.*).tsv'),
		 r'\1/collapsed/\2-\3.tsv')

def collapseIsoformRbpBinding(infiles, outfile):

	# Run
	run_r_job('collapse_isoform_rbp_binding', infiles, outfile, modules=['R/4.0.3'], W='00:30', GB=30, n=1, run_locally=True, print_outfile=False, print_cmd=False)

#############################################
########## 6. Get BED file
#############################################

@collate('arion/illumina/s11-motif_scan.dir/human/intersect/*/*-isoform_intersect_coordinates.bed',
		 regex(r'(.*)/(.*)/intersect.*/.*/(.*?_.*?)_.*-(.*).bed'),
		 r'\1/\2/rbp_intersect_bed/\3-cluster_isoform_overlap.bed')

def getIsoformRbpBindingBed(infiles, outfile):
	
	# Command
	infiles_str = ' '.join(infiles)
	cmd_str = ''' cat {infiles_str} > {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, run_locally=True, print_outfile=False)


#############################################
########## 6. Get AS exon bed
#############################################

@collate('arion/illumina/s08-suppa.dir/human/05-summaries/human-*-suppa_summary.tsv',
		   regex(r'(.*)/s08-suppa.dir/(.*?)/.*/(.*?)-.*-suppa_summary.tsv'),
		   r'\1/s11-motif_scan.dir/human/bed/\2_SE_500bp.bed')

def getAsExonBed(infiles, outfile):

	# Run
	run_r_job('get_as_exon_bed', infiles, outfile, modules=['R/4.0.3'], W='00:15', GB=10, n=1, run_locally=True, print_outfile=False, print_cmd=False)

#############################################
########## 7. Intersect AS exons and motifs
#############################################

@transform('arion/illumina/s11-motif_scan.dir/human/homer/*/*_homer_results.txt',
		   regex(r'(.*)/homer/(.*)/.*.txt'),
		   add_inputs(getAsExonBed),
		   r'\1/intersect/\2/\2-SE_intersect.tsv')

def intersectAsExonMotifs(infiles, outfile):

	# Run
	run_r_job('intersect_as_exon_motifs', infiles, outfile, modules=['R/4.0.3'], W='00:30', GB=30, n=1, run_locally=False, print_outfile=True, print_cmd=False)

#############################################
########## 8. Convert to BED
#############################################

@collate('arion/illumina/s11-motif_scan.dir/human/intersect/*/*-SE_intersect.tsv',
		 regex(r'(.*)/intersect/.*?/.*-(.*_intersect).tsv'),
		 r'\1/bed/\2.bed')

def mergeAsExonMotifs(infiles, outfile):

	# Run
	run_r_job('convert_rbp_scan_bed', infile, outfile, modules=['R/4.0.3'], W='00:15', GB=10, n=1, run_locally=True, print_outfile=False, print_cmd=False)

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