#################################################################
#################################################################
############### GEO Illumina RNA-Seq ################
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
import pandas as pd
import numpy as np
# from rpy2.robjects import r, pandas2ri
# pandas2ri.activate()

##### 2. LSF #####
# 2.1 Import
sys.path.append('/hpc/users/torred23/pipelines/support')
import lsf

# 2.2 Default parameters
r_source = 'pipeline/scripts/geo-illumina.R'
py_source = 'pipeline/scripts/GeoIllumina.py'
P = 'acc_apollo'
q = 'premium'
W = '00:30'
GB = 5
n = 1
mkdir_val = True

# 2.3 Wrappers
# CMD
def run_job(cmd_str, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_job(cmd_str, outfile, P=P, W = W, GB = GB, n = n, q = q, mkdir=mkdir_val, **kwargs)

# R
def run_r_job(func_name, func_input, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_r_job(func_name, func_input, outfile, r_source=r_source, P=P, W = W, GB = GB, n = n, q = q, mkdir=mkdir_val, **kwargs)

# Py
def run_py_job(func_name, func_input, outfile, W = W, GB = GB, n = n, q = q, **kwargs):
	lsf.run_py_job(func_name, func_input, outfile, py_source=py_source, P=P, W = W, GB = GB, n = n, q = q, mkdir=mkdir_val, **kwargs)

##### 3. Custom script imports #####
# 3.1 Python
#sys.path.append('pipeline/scripts')
#import GeoIllumina as P

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
with open('/sc/arion/projects/GuccioneLab/genome-indices/genome-indices.json') as openfile:
	genome_indices = json.load(openfile)

##### 3. Variables #####
linked_fastq = glob.glob('arion/geo_illumina/s01-datasets.dir/*/fastq/*/*.fastq.gz')
human_linked_fastq = [x for x in linked_fastq if 'wang' not in x and 'boroviak' not in x]
primate_linked_fastq = {
	'macaque': [x for x in linked_fastq if 'wang' in x],
	'marmoset': [x for x in linked_fastq if 'boroviak' in x]
}
human_filtered_gtf = 'arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.gtf'
human_star_index = 'arion/illumina/s03-alignment.dir/human/isoseq/STAR/index'
human_rsem_index = 'arion/illumina/s03-alignment.dir/human/isoseq/RSEM/index/Homo_sapiens.GRCh38.102_talon-SJ_filtered.idx.fa'#arion/illumina/s04-alignment.dir/human/all/RSEM/index/Homo_sapiens.GRCh38.102_talon-all-SJ_filtered.idx.fa
human_junction_file = 'to_fix' #arion/isoseq/s05-talon.dir/human/Homo_sapiens.GRCh38.102_talon_junctions.tsv
nonhuman_genomes = {
	'macaque': {
		'gtf': 'arion/datasets/reference_genomes/macaque/Macaca_mulatta.Mmul_10.102.gtf',
		'gtf_lifted': 'arion/geo_illumina/s05-primates.dir/macaque/gtf/macaque-hg38_filtered_lifted.gtf',
		'genome_fasta': 'arion/datasets/reference_genomes/macaque/Macaca_mulatta.Mmul_10.dna_sm.primary_assembly.fa'#,
		# 'sashimi_bams': 'arion/geo_illumina/summary.dir/sashimi/settings/macaque-bams.txt'
	},
	'marmoset': {
		'gtf': 'arion/datasets/reference_genomes/marmoset/ncbiRefSeq.gtf',
		'gtf_lifted': 'arion/geo_illumina/s05-primates.dir/marmoset/gtf/marmoset-hg38_filtered_lifted.gtf',
		'genome_fasta': 'arion/datasets/reference_genomes/marmoset/calJac4.fa'#,
		# 'sashimi_bams': 'arion/geo_illumina/summary.dir/sashimi/settings/marmoset-bams_subset.txt'
	}
}
# macaque_gtf = 'arion/datasets/reference_genomes/macaque/Macaca_mulatta.Mmul_10.102.gtf'
# macaque_primary_assembly = 'arion/datasets/reference_genomes/macaque/Macaca_mulatta.Mmul_10.dna_sm.primary_assembly.fa'
# marmoset_gtf = 'arion/datasets/reference_genomes/marmoset/ncbiRefSeq.gtf'
# marmoset_primary_assembly = 'arion/datasets/reference_genomes/marmoset/calJac4.fa'

#######################################################
#######################################################
########## S1. Fix sample names
#######################################################
#######################################################

#############################################
########## 1. Yan
#############################################

@transform('arion/datasets/yan/yan-samples.csv',
		   regex(r'.*/(.*)/(.*)s.csv'),
		   r'arion/geo_illumina/s01-datasets.dir/\1/\2_names.csv')

def fixYanSamples(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
	
	# Columns to rename
	column_dict = {'Run': 'Run', 'Assay Type': 'assay_type', 'GEO_Accession (exp)': 'geo_sample_id', 'Cell_type': 'geo_cell_type', 'Development_stage': 'geo_development_stage', 'Library Name': 'geo_library_name'}

	# Read
	sample_dataframe = pd.read_csv(infile, comment='#')

	# Rename and select embryo samples
	sample_dataframe = sample_dataframe.rename(columns=column_dict)[column_dict.values()].query('geo_development_stage == "early blastomere"')

	# Add cell type
	sample_dataframe['cell_type'] = [x.split(':')[-1].replace('-cell embryo', 'C').replace('Morulae', 'morula').replace('Late blastocyst', 'blastocyst').replace('Zygote', '1C').replace('-cell', 'C').replace('Oocyte', 'oocyte').replace(' ', '') for x in sample_dataframe['geo_library_name']]

	# Sort and add sample number
	sample_dataframe = sample_dataframe.sort_values('cell_type')
	sample_dataframe['sample_number'] = sample_dataframe.groupby('cell_type').cumcount()+1

	# Add sample name
	sample_dataframe['sample_name'] = ['human_{cell_type}_Rep{sample_number}'.format(**rowData) for index, rowData in sample_dataframe.iterrows()]

	# Write
	sample_dataframe.to_csv(outfile, index=False)

#############################################
########## 2. Yan
#############################################

@transform('arion/datasets/xue/xue-samples.csv',
		   regex(r'.*/(.*)/(.*)s.csv'),
		   r'arion/geo_illumina/s01-datasets.dir/\1/\2_names.csv')

def fixXueSamples(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)
		
	# Columns to rename
	column_dict = {'Run': 'Run', 'Assay Type': 'assay_type', 'GEO_Accession (exp)': 'geo_sample_id', 'Cell_type': 'geo_cell_type', 'Organism': 'geo_organism', 'source_name': 'geo_source_name'}

	# Read
	sample_dataframe = pd.read_csv(infile, comment='#')

	# Rename and select embryo samples
	sample_dataframe = sample_dataframe.rename(columns=column_dict)[column_dict.values()].query('geo_organism == "Homo sapiens" and geo_cell_type != "primary blood"')

	# Add cell type
	sample_dataframe['cell_type'] = [x.split(':')[-1].replace('-cell blastomere', 'C').replace('zygote', '1C') for x in sample_dataframe['geo_cell_type']]
	sample_dataframe

	# Sort and add sample number
	sample_dataframe = sample_dataframe.sort_values('cell_type')
	sample_dataframe['sample_number'] = sample_dataframe.groupby('cell_type').cumcount()+1

	# Add sample name
	sample_dataframe['sample_name'] = ['human_{cell_type}_Rep{sample_number}'.format(**rowData) for index, rowData in sample_dataframe.iterrows()]

	# Write
	sample_dataframe.to_csv(outfile, index=False)

#############################################
########## 3. Liu
#############################################

@transform('arion/datasets/liu/liu-samples.csv',
		   regex(r'.*/(.*)/(.*)s.csv'),
		   r'arion/geo_illumina/s01-datasets.dir/\1/\2_names.csv')

def fixLiuSamples(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Columns to rename
	column_dict = {'Run': 'Run', 'Assay Type': 'assay_type', 'Sample Name': 'geo_sample_name', 'Organism': 'geo_organism', 'Library Name': 'geo_library_name'}

	# Read
	sample_dataframe = pd.read_csv(infile, comment='#')

	# Rename and select embryo samples
	sample_dataframe = sample_dataframe.rename(columns=column_dict)[column_dict.values()].query('geo_organism == "Homo sapiens" and assay_type == "RNA-Seq"')

	# Add sample name
	sample_dataframe['sample_name'] = ['human_'+x.replace('R_', '').replace('-cell', 'C').replace('Morula', 'morula').replace('Zygote', '1C') for x in sample_dataframe['geo_sample_name']]

	# Write
	sample_dataframe.to_csv(outfile, index=False)

#############################################
########## 4. Wang
#############################################

@transform('arion/datasets/wang/wang-samples.csv',
		   regex(r'.*/(.*)/(.*)s.csv'),
		   r'arion/geo_illumina/s01-datasets.dir/\1/\2_names.csv')

def fixWangSamples(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Columns to rename
	column_dict = {'Run': 'Run', 'source_name': 'source_name'}

	# Read
	sample_dataframe = pd.read_csv(infile, comment='#')

	# Rename and select embryo samples
	sample_dataframe = sample_dataframe.rename(columns=column_dict)[column_dict.values()]

	# Add sample name
	sample_dict = {'8-cell embryo': '8C', 'blastocyst': 'blastocyst', 'morula': 'morula', 'germinal vesicle-stage oocyte': 'GV_oocyte', '1-cell embryo at pronuclear (PN) stage': '1C', 'mature oocyte': 'M_oocyte', '2-cell embryo': '2C', '4-cell embryo': '4C'}
	sample_dataframe['cell_type'] = [sample_dict[x] for x in sample_dataframe['source_name']]

	# Sort and add sample number
	sample_dataframe = sample_dataframe.sort_values('cell_type')
	sample_dataframe['sample_number'] = sample_dataframe.groupby('cell_type').cumcount()+1

	# Add sample name
	sample_dataframe['sample_name'] = ['macaque_{cell_type}_Rep{sample_number}'.format(**rowData) for index, rowData in sample_dataframe.iterrows()]

	# Write
	sample_dataframe.to_csv(outfile, index=False)

#############################################
########## 5. Boroviak
#############################################

@transform('arion/datasets/boroviak/boroviak-samples.tsv',
		   regex(r'.*/(.*)/(.*)s.tsv'),
		   r'arion/geo_illumina/s01-datasets.dir/\1/\2_names.csv')

def fixBoroviakSamples(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Columns to rename
	column_dict = {'run_accession': 'Run', 'sample_title': 'sample_title'}

	# Read
	sample_dataframe = pd.read_table(infile, comment='#')

	# Rename and select embryo samples
	sample_dataframe = sample_dataframe.rename(columns=column_dict)[column_dict.values()]

	# Add sample name
	sample_dataframe['sample_name'] = ['marmoset_'+x.replace('Zygote', '1C').replace('cell', 'C') for x in sample_dataframe['sample_title']]

	# Write
	sample_dataframe.to_csv(outfile, index=False)

#############################################
########## 6. Link
#############################################

@subdivide((fixYanSamples, fixXueSamples, fixLiuSamples, fixWangSamples, fixBoroviakSamples), #fixBoroviakSamples
		   regex(r'(.*)/.*-sample_names.csv'),
		   r'\1/fastq/*/*.fastq.gz',
		   r'\1/fastq/{sample_name}/{fastq_basename}')

def linkFASTQ(infile, outfiles, outfileRoot):

	# Read
	sample_dataframe = pd.read_csv(infile)

	# Get dataset name
	dataset = os.path.basename(infile).split('-')[0]

	# Loop through dataframe
	for index, rowData in sample_dataframe.iterrows():

		# Get FASTQ files
		fastq_files = glob.glob('/hpc/users/torred23/pipelines/projects/early-embryo/geo-illumina/arion/datasets/{dataset}/rawdata/{Run}*.fastq.gz'.format(**locals(), **rowData))

		# Loop
		for fastq_file in fastq_files:

			# Get basename
			fastq_basename = os.path.basename(fastq_file)

			# Get outfile
			outfile = outfileRoot.format(**locals(), **rowData)
		
			# Outdir
			outdir = os.path.dirname(outfile)
			if not os.path.exists(outdir):
				os.makedirs(outdir)

			# Link
			if not os.path.exists(outfile):
				os.system('ln -s {fastq_file} {outfile}'.format(**locals()))

#######################################################
#######################################################
########## S2. FASTQC
#######################################################
#######################################################

#############################################
########## 1. Run
#############################################

# @follows(linkFASTQ)

@transform((linked_fastq),
		   regex(r'(.*)/s01-datasets.dir/(.*)/fastq/(.*)/(.*).fastq.gz'),
		   r'\1/s02-fastqc.dir/\2/\4_fastqc.html')

def runFastQC(infile, outfile):

	# Command
	cmd_str = '''fastqc --outdir=$(dirname {outfile}) {infile}'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['fastqc/0.11.8'], W='03:00', GB=12, n=1, print_outfile=False)

#############################################
########## 2. MultiQC
#############################################

# @follows(runFastQC)

@transform(('arion/geo_illumina/s02-fastqc.dir/*', 'arion/geo_illumina/s03-alignment.dir/*', 'arion/geo_illumina/s05-primates.dir/*'),
		   regex(r'(.*)/s..-(.*).dir/(.*)'),
		   r'\1/multiqc/\3/\3_\2/multiqc_report.html')

def runMultiQC(infile, outfile):

	# Command
	cmd_str = 'multiqc --outdir $(dirname {outfile}) {infile}'.format(**locals())

	# Run
	if not os.path.exists(outfile):
		run_job(cmd_str, outfile, conda_env='env', W="01:00", GB=10, n=1, print_outfile=False, run_locally=False, stdout=outfile.replace('.html', '.log'))

#######################################################
#######################################################
########## S3. Align
#######################################################
#######################################################

#############################################
########## 1. STAR junctions
#############################################

# @follows(linkFASTQ)

# @collate(human_linked_fastq, # 3 datasets
# 		 regex(r'(.*)/s01-datasets.dir/(.*)/fastq/(.*)/.*.fastq.gz'),
# 		 add_inputs(human_star_index),
# 		 r'\1/s03-alignment.dir/\2/STAR/pass1/\3/\3-SJ.out.tab')

@collate('arion/datasets/asami/rawdata/*.fastq.gz', # asami
		 regex(r'arion/datasets/(.*)/rawdata/(.*?)_..fastq.gz'),
		 add_inputs(human_star_index),
		 r'arion/geo_illumina/s03-alignment.dir/\1/STAR/pass1/\2/\2-SJ.out.tab')

def getStarJunctions(infiles, outfile):

	# Split
	fastq_files = [x[0] for x in infiles]
	star_index = infiles[0][1]

	# FASTQ string
	fastq_str = ' '.join(fastq_files)

	# Prefix
	prefix = outfile[:-len('SJ.out.tab')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_str} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 100 \
		--outSAMtype None'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="06:00", print_cmd=True, GB=5, n=15, modules=['star/2.7.5b'], stdout=outfile.replace('-SJ.out.tab', '_job.log'))

# ls arion/geo_illumina/s03-alignment.dir/*/STAR/pass1/*/*_job.log | jsc

#############################################
########## 2. STAR BAM
#############################################

# @follows(getStarJunctions)

# @collate(human_linked_fastq, # 3 datasets
# 		 regex(r'(.*)/s01-datasets.dir/(.*)/fastq/(.*)/.*.fastq.gz'),
# 		 add_inputs(human_star_index, r'\1/s03-alignment.dir/\2/STAR/pass1/*/*-SJ.out.tab'),
# 		 r'\1/s03-alignment.dir/\2/STAR/pass2/\3/\3-Aligned.sortedByCoord.out.bam')

@collate('arion/datasets/asami/rawdata/*.fastq.gz', # asami
		 regex(r'arion/datasets/(.*)/rawdata/(.*?)_..fastq.gz'),
		 add_inputs(human_star_index, r'arion/geo_illumina/s03-alignment.dir/\1/STAR/pass1/*/*-SJ.out.tab'),
		 r'arion/geo_illumina/s03-alignment.dir/\1/STAR/pass2/\2/\2-Aligned.sortedByCoord.out.bam')

def runStar(infiles, outfile):

	# Split
	fastq_files = [x[0] for x in infiles]
	star_index = infiles[0][1]
	sj_files = infiles[0][2:]

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
	fastq_str = ' '.join(fastq_files)
	sj_files_str = ' '.join(sj_files)

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_str} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 50000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="06:00", GB=15, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# ls arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*.log | jsc
# ls arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*Log.out | lr
# ls arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*.err | xargs wc -l

# ls arion/geo_illumina/s03-alignment.dir/asami/STAR/pass2/*/*.log | jsc
# ls arion/geo_illumina/s03-alignment.dir/asami/STAR/pass2/*/*Log.out | lr
# ls arion/geo_illumina/s03-alignment.dir/asami/STAR/pass2/*/*.err | xargs wc -l

#############################################
########## 3. Junction counts
#############################################

# @follows(runStar)

# @collate('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*-SJ.out.tab',
# 		 regex(r'(.*)/(.*)/(.*)/pass2/.*/.*.tab'),
# 		 add_inputs(human_junction_file),
# 		 r'\1/\2/\3/\2-junction_counts.tsv')

# def getJunctionCounts(infiles, outfile):

# 	# Infiles
# 	sj_files = [x[0] for x in infiles]
# 	junction_file = infiles[0][1]

# 	# Run
# 	run_r_job('get_junction_counts', sj_files, outfile, additional_params=junction_file, print_outfile=False, conda_env='env', W='02:00', GB=50, n=3, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

# find arion/geo_illumina/s03-alignment.dir -name "*junction_counts*" | xargs rm

#############################################
########## 4. RSEM expression
#############################################

# @follows(runStar)

# # @transform('arion/illumina/s04-alignment.dir/human/2C_vs_4C/STAR/pass2/human_4C_B3_9/human_4C_B3_9-Aligned.toTranscriptome.out.bam',
# @transform('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/human_morula_Rep15-Aligned.toTranscriptome.out.bam',
@transform('arion/geo_illumina/s03-alignment.dir/xue/STAR/pass2/*/*-Aligned.toTranscriptome.out.bam',
# @transform('arion/geo_illumina/s03-alignment.dir/xue/STAR/pass2/human_1C_Rep1/*-Aligned.toTranscriptome.out.bam',
		   regex(r'(.*)/STAR/.*/(.*)-Aligned.toTranscriptome.out.bam'),
		   add_inputs(human_rsem_index),
		   r'\1/RSEM_forward/\2/\2.isoforms.results')

def runRsem(infiles, outfile):

	# Variables
	prefix = outfile[:-len('.isoforms.results')]
	reference_name = infiles[1][:-len('.idx.fa')]

	# Paired end
	dataset = outfile.split('/')[-4]
	paired_str = '--paired-end' if dataset in ['xue', 'asami'] else ''
	strandedness_str = 'forward'#'forward' if dataset == 'xue' else 'none'

	# Command
		# --strandedness reverse \
	cmd_str = '''rsem-calculate-expression \
		--alignments \
		--strandedness {strandedness_str} \
		{paired_str} \
		--estimate-rspd \
		--num-threads 200 \
		--no-bam-output \
		{infiles[0]} \
		{reference_name} \
		{prefix} > {prefix}.rsem.log && \
		rsem-plot-model {prefix} {prefix}.quant.pdf '''.format(**locals())
		# --calc-ci \

	# Run
	run_job(cmd_str, outfile, W="06:00", GB=2, n=25, modules=['rsem/1.3.3'], print_cmd=False, stdout=outfile.replace('.isoforms.results', '.log'), stderr=outfile.replace('.isoforms.results', '.err'))

# find arion/geo_illumina/s03-alignment.dir/*/RSEM -name "*.log" | grep -v 'rsem' | jsc

# ml salmon/1.2.1
# salmon quant -i /sc/arion/projects/GuccioneLab/genome-indices/hg38/salmon/index -l A -1 /hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/asami/rawdata/SRR12624495_1.fastq.gz -2 /hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/asami/rawdata/SRR12624495_2.fastq.gz -p 8 --validateMappings -o /sc/arion/scratch/torred23/salmon/SRR12624495/ &> SRR12624495.salmon.log
# asami ISF?

#############################################
########## 4. Merge RSEM expression
#############################################

@transform('arion/geo_illumina/s07-scrna.dir/RSEM_individual/*',
		   regex(r'(.*)/RSEM_individual/(.*)'),
		   r'\1/RSEM_grouped/\2_average_grouped_expression.tsv')

def mergeRsem(infile, outfile):

	# Run
	run_r_job('merge_rsem', infile, outfile, modules=['R/4.0.3'], W='03:00', GB=50, n=5, run_locally=False, print_outfile=False, print_cmd=False, q='sla')

# rm -r /hpc/users/torred23/pipelines/projects/early-embryo/arion/geo_illumina/s07-scrna.dir/RSEM_grouped

#############################################
########## 6. Create BigWig
#############################################

# @transform(runStar,
# 		   suffix('-Aligned.sortedByCoord.out.bam'),
# 		   '.bw')

# def createBigWig(infile, outfile):

# 	# Command
# 	cmd_str = """bamCoverage --outFileFormat=bigwig --skipNonCoveredRegions --numberOfProcessors=50 --normalizeUsing RPKM -b {infile} -o {outfile}""".format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', W='05:00', n=7, GB=6)

#######################################################
#######################################################
########## S5. Expression
#######################################################
#######################################################

#############################################
########## 1. Prepare metadata
#############################################

# @collate('arion/geo_illumina/s03-alignment.dir/*/RSEM/*/*.isoforms.results',
# @collate(runRsem,
# 		 regex(r'(.*)/s03-alignment.dir/(.*)/RSEM/.*.isoforms.results'),
# 		 r'\1/s04-expression.dir/\2/\2-sample_metadata.txt')
@collate('arion/geo_illumina/s03-alignment.dir/xue/RSEM_*/*/*.isoforms.results',
		 regex(r'(.*)/s03-alignment.dir/(.*)/RSEM_(.*?)/.*.isoforms.results'),
		 r'\1/s04-expression.dir/\2_\3/\2_\3-sample_metadata.txt')

def prepareSampleMetadata(infiles, outfile):

	# Create directory
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Prepare dict
	sample_dict = [{'files': x, 'names': x.split('/')[-2], 'cell_type': x.split('/')[-2].split('_')[1]} for x in infiles]

	# Convert to dataframe
	sample_dataframe = pd.DataFrame(sample_dict)

	# Write
	sample_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 2. Aggregate
#############################################

@transform(prepareSampleMetadata,
		   suffix('-sample_metadata.txt'),
		   add_inputs(human_filtered_gtf),
		   '-counts.rda')

def aggregateCounts(infiles, outfile):

	# Run
	run_r_job('aggregate_counts', infiles, outfile, conda_env='env', run_locally=False, W='00:30', GB=30, n=1, wait=False, ow=True, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 3. Get size factors
#############################################

@transform(('arion/geo_illumina/s04-expression.dir/*/*-counts.rda', 'arion/geo_illumina/s05-primates.dir/*/RSEM/counts/*-counts.rda'),
		   suffix('-counts.rda'),
		   '-size_factors.tsv')

def getSizeFactors(infile, outfile):

	# Run
	run_r_job('get_size_factors', infile, outfile, conda_env='env', modules=[], W='00:15', GB=20, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 4. Create scaled BigWig
#############################################

@transform('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
		   regex(r'(.*.dir)/(.*?)/(.*)-Aligned.sortedByCoord.out.bam'),
		   add_inputs(r'arion/geo_illumina/s04-expression.dir/\2/\2-size_factors.tsv'),
		   r'\1/\2/\3-scaled_3.bw')

def createScaledBigWig(infiles, outfile):

	# Read size factor
	normalization_dict = pd.read_table(infiles[1], index_col='sample_name')['size_factor_reciprocal'].to_dict()
	size_factor = normalization_dict[outfile.split('/')[-2]]

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=3 --skipNonCoveredRegions --numberOfProcessors=48 --scaleFactor {size_factor} -b {infiles[0]} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='02:00', n=8, GB=4, print_outfile=False)

# find arion/geo_illumina/s03-alignment.dir/liu/STAR/pass2 -name "*scaled.bw"

#############################################
########## 7. Merge scaled BigWig
#############################################

# @collate(createScaledBigWig,
@collate('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*-scaled_3.bw',
		 regex(r'(.*)/s03-alignment.dir/(.*)/STAR/pass2/(.*?)_Rep.*/*-scaled_3.bw'),
		 add_inputs(r'arion/datasets/reference_genomes/human/*.chrom.sizes'),
		 r'\1/s04-expression.dir/\2/scaled_bw/\3-\2.bw')

def mergeScaledBigWig(infiles, outfile):
	
	# Files
	wig_file = outfile.replace('.bw', '.wig')
	bedgraph_file = outfile.replace('.bw', '.bedgraph')
	infiles_str = ' '.join([x[0] for x in infiles])

	# Command
	# cmd_str = """ wiggletools mean {infiles_str} > {wig_file} && wigToBigWig {wig_file} {infiles[0][1]} {outfile} && rm {wig_file} """.format(**locals())
	cmd_str = """ wiggletools mean {infiles_str} | sed '/^KI.*/{{s///;q;}}' > {wig_file} && wigToBigWig {wig_file} {infiles[0][1]} {outfile} && rm {wig_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} > {wig_file} && wiggletools write_bg - {wig_file} | sort -k1,1 -k2,2n > {bedgraph_file} && bedGraphToBigWig {bedgraph_file} {infiles[0][1]} {outfile} && rm {wig_file} {bedgraph_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} | sed '/^KI.*/{{s///;q;}}' > {wig_file} && wiggletools write_bg - {wig_file} | sort -k1,1 -k2,2n > {bedgraph_file} && bedGraphToBigWig {bedgraph_file} {infiles[0][1]} {outfile} && rm {wig_file} {bedgraph_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} > {wig_file} && wiggletools write_bg - {wig_file} | sort -k1,1 -k2,2n > {bedgraph_file} && bedGraphToBigWig {bedgraph_file} {infiles[0][1]} {outfile} && rm {wig_file} {bedgraph_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} > {wig_file} """.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['wiggletools/1.2', 'ucsc-utils/2020-03-17'], W='00:30', n=1, GB=10, print_outfile=False, stdout=outfile.replace('.bw', '.log'), stderr=outfile.replace('.bw', '.err'))

#######################################################
#######################################################
########## S5. Sashimi
#######################################################
#######################################################

#############################################
########## 1. BAM tables
#############################################

# @collate(runStar,
# @collate('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
# @collate('arion/illumina/s04-alignment.dir/human/all/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
@collate(('arion/geo_illumina/s03-alignment.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam', 'arion/illumina/s04-alignment.dir/human/all/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam', 'arion/geo_illumina/s05-primates.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam'),
		 regex(r'.*.dir/(.*)/STAR/pass2/.*/.*.bam'),
		 r'arion/geo_illumina/summary.dir/sashimi/settings/\1-bams.txt')

def makeBamTables(infiles, outfile):
		
	# Create dataframe
	bam_dataframe = pd.DataFrame([{'sample_name': x.split('/')[-2], 'bam': os.path.join(os.getcwd(), x), 'cell_type': x.split('/')[-2].split('_')[1]} for x in infiles])

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Write
	bam_dataframe.to_csv(outfile, sep='\t', index=False, header=False)

#############################################
########## 2. Sashimi plots
#############################################

# def sashimiJobs():

# 	# Loop through comparisons
# 	for comparison, comparison_info in comparison_dict.items():
		
# 		# Get BAM paths
# 		bams = {}
# 		for group in comparison_info['contrast']:
# 			with open('arion/illumina/s07-rmats.dir/isoseq/{comparison}/bams/{group}-bams.txt'.format(**locals())) as openfile:
# 				bams[group] = [os.path.join(os.getcwd(), x) for x in openfile.read().split(',')]
				
# 		# Read splicing results
# 		splicing_dataframe = pd.read_table('arion/illumina/summary.dir/{comparison}/{comparison}-isoseq-splicing_summary.tsv'.format(**locals()))
# 		splicing_dataframe['i'] = splicing_dataframe['pval_rmats'].sort_values().index

# 		# Filter
# 		splicing_dataframe = splicing_dataframe.query('significant_rmats==True and significant_suppa==True and event_type == "SE"').sort_values('pval_rmats')
		
# 		# Loop through splicing results
# 		for index, rowData in splicing_dataframe.iterrows():
# 			min_coverage = 5
# 			plot_settings = {
# 				'groups': [{'name': group, 'bams': bams[group], 'psi': round(rowData['PSI'+str(i+1)], ndigits=2), 'color': '#377eb8' if i == 0 else '#e41a1c'} for i, group in enumerate(comparison_info['contrast'])],
# 				'event_id': rowData['Event_id'],
# 				'gene_name': rowData['gene_name'],
# 				'rmats_statistics': 'dPSI={dPSI_rmats:.2f}, FDR={pval_rmats:.1e}'.format(**rowData),
# 				'gtf': reference_dict['isoseq']['gtf_cds'],
# 				'min_coverage': min_coverage
# 			}

# 			# Files
# 			infile = 'arion/illumina/summary.dir/{comparison}/plots/sashimi/{gene_name}-{event_type}-{i}-{min_coverage}.json'.format(**locals(), **rowData)
# 			outfile = infile.replace('.json', '.pdf')
			
# 			# Outdir
# 			outdir = os.path.dirname(outfile)
# 			if not os.path.exists(outdir):
# 				os.makedirs(outdir)

# 			# Write infile
# 			with open(infile, 'w') as openfile:
# 				json.dump(plot_settings, openfile)

# 			yield [infile, outfile]


# # @follows(functionToFollow)

# @files(sashimiJobs)

# def sashimiPlot(infile, outfile):

# 	# Run
# 	run_py_job('ggsashimi', infile, outfile, W="00:10", GB=30, n=1, run_locally=False, conda_env='env')

#############################################
########## 3. Novel primate genes
#############################################

@collate('arion/geo_illumina/s05-primates.dir-old/*/RSEM/counts/*-counts.rda',
		 regex(r'(.*)/s05.*.rda'),
		 r'\1/summary.dir/sashimi-novel/primate-novel_genes.tsv')

def getNovelPrimateGenes(infiles, outfile):

	# Run
	run_r_job('get_novel_primate_genes', infiles, outfile, conda_env='env')

#############################################
########## 4. Novel primate sashimi
#############################################

def novelGeneJobs():
	gene_dataframe = pd.read_table('arion/geo_illumina/summary.dir/sashimi-novel/primate-novel_genes.tsv')
	gene_ids = gene_dataframe.query('macaque > 3 and marmoset > 3')['gene_id']#[:10]
	sashimi_genomes = nonhuman_genomes.copy()
	sashimi_genomes['human'] = {'gtf': 'arion/illumina/s04-alignment.dir/human/all/gtf/Homo_sapiens.GRCh38.102_talon-all-SJ_filtered.gtf', 'sashimi_bams':'arion/geo_illumina/summary.dir/sashimi/settings/all-bams.txt'}
	for gene_id in gene_ids:
		for organism, parameter_dict in sashimi_genomes.items():
			if 'gtf_lifted' in parameter_dict.keys():
				parameter_dict['gtf'] = parameter_dict['gtf_lifted']
			infiles = [parameter_dict['gtf'], parameter_dict['sashimi_bams']]
			outfile = 'arion/geo_illumina/summary.dir/sashimi-novel/plots/split/{organism}-{gene_id}.png'.format(**locals())
			yield [infiles, outfile, gene_id]

# @follows(getNovelPrimateGenes)

@files(novelGeneJobs)

def novelPrimateGeneSashimi(infiles, outfile, gene_id):

	# Plot parameters
	plot_prefix, plot_format = outfile.rsplit('.', 1)

	# Command
	cmd_str = ''' /sc/arion/work/torred23/libraries/ggsashimi/ggsashimi-gene_id/ggsashimi.py \
		--bam {infiles[1]} \
		--gtf {infiles[0]} \
		--gene_id {gene_id} \
		--min-coverage 1000 \
		--overlay 3 \
		--fix-y-scale \
		--aggr mean_j \
		--out-prefix {plot_prefix} \
		--out-format {plot_format}
	'''.format(**locals())

	# Run
	# if 'marmoset' in outfile:
	run_job(cmd_str, outfile, conda_env='env', modules=['samtools/1.9'], W='03:00', GB=25, n=1, print_cmd=False)#, print_cmd=False)

#######################################################
#######################################################
########## S5. Primate analysis
#######################################################
#######################################################

#############################################
########## 1. Copy lifted GTF
#############################################

def gtfJobs():
	gtf_dict = {
		'macaque': 'arion/isoseq/s10-liftover.dir/human/merged/hg38ToRheMac10/Homo_sapiens.GRCh38.102_talon-SJ_filtered.cds-hg38ToRheMac10-gene_id.gtf',
		'marmoset': 'arion/isoseq/s10-liftover.dir/human/merged/hg38ToCalJac4/Homo_sapiens.GRCh38.102_talon-SJ_filtered.cds-hg38ToCalJac4-gene_id.gtf'
	}
	for organism, infile in gtf_dict.items():
		outfile = 'arion/geo_illumina/s05-primates.dir/{organism}/gtf/{organism}-hg38_filtered_lifted.gtf'.format(**locals())
		yield [infile, outfile]

@files(gtfJobs)

def copyLiftedGTF(infile, outfile):

	# Run
	run_r_job('copy_lifted_gtf', infile, outfile, conda_env='env', W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 2. STAR index
#############################################

@transform(copyLiftedGTF,
		   regex(r'(.*)/(.*)/gtf/.*.gtf'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.fa'),
		   r'\1/\2/STAR/index')

def buildStarIndex(infiles, outfile):

	# Filter infiles
	infiles = [x for x in infiles if 'cdna' not in x and 'renamed' not in x]

	# Command
	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {infiles[1]} --sjdbGTFfile {infiles[0]} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='02:00', GB=5, n=15, ow=True, print_cmd=True, stdout=os.path.join(outfile, 'job.log'), jobname='_'.join(outfile.split('/')[-3:]), wait=False)

#############################################
########## 3. RSEM index
#############################################

@transform(copyLiftedGTF,
		   regex(r'(.*)/(.*)/gtf/.*.gtf'),
		   add_inputs(r'arion/datasets/reference_genomes/\2/*.fa'),
		   r'\1/\2/RSEM/index/\2-hg38_filtered_lifted.idx.fa')

def buildRsemIndex(infiles, outfile):

	# Filter infiles
	infiles = [x for x in infiles if 'cdna' not in x and 'renamed' not in x]

	# Command
	basename = outfile[:-len('.idx.fa')]
	cmd_str = ''' rsem-prepare-reference --gtf {infiles[0]} --num-threads 10 {infiles[1]} {basename} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="00:30", GB=5, n=3, modules=['rsem/1.3.3'], print_cmd=False, stdout=basename+'.log', stderr=basename+'.err')

#############################################
########## 4. STAR first pass
#############################################

def primateSjJobs():
	for organism, fastq_files in primate_linked_fastq.items():
		fastq_dataframe = pd.DataFrame({'fastq_path': fastq_files}).sort_values('fastq_path')
		fastq_dataframe['sample_name'] = [x.split('/')[-2] for x in fastq_dataframe['fastq_path']]
		fastq_dict = fastq_dataframe.groupby('sample_name')['fastq_path'].apply(lambda x: list(x)).to_dict()
		for sample_name, fastq_pair in fastq_dict.items():
			infiles = fastq_pair+['arion/geo_illumina/s05-primates.dir/{organism}/STAR/index'.format(**locals())]
			outfile = 'arion/geo_illumina/s05-primates.dir/{organism}/STAR/pass1/{sample_name}/{sample_name}-SJ.out.tab'.format(**locals())
			yield [infiles, outfile]

@follows(buildStarIndex)

@files(primateSjJobs)

def getPrimateStarJunctions(infiles, outfile):

	# Split
	fastq_files = infiles[:-1]
	star_index = infiles[-1]

	# FASTQ string
	fastq_str = ' '.join(fastq_files)

	# Prefix
	prefix = outfile[:-len('SJ.out.tab')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_str} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 100 \
		--outSAMtype None'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="02:00", print_outfile=True, GB=5, n=15, modules=['star/2.7.5b'], stdout=outfile.replace('-SJ.out.tab', '_job.log'))

# find arion/geo_illumina/s05-primates.dir/*/STAR/pass1/ -name "*job.log" | jsc

#############################################
########## 5. STAR second pass
#############################################

def primateStarJobs():
	for organism, fastq_files in primate_linked_fastq.items():
		fastq_dataframe = pd.DataFrame({'fastq_path': fastq_files}).sort_values('fastq_path')
		fastq_dataframe['sample_name'] = [x.split('/')[-2] for x in fastq_dataframe['fastq_path']]
		fastq_dict = fastq_dataframe.groupby('sample_name')['fastq_path'].apply(lambda x: list(x)).to_dict()
		for sample_name, fastq_pair in fastq_dict.items():
			star_index = 'arion/geo_illumina/s05-primates.dir/{organism}/STAR/index'.format(**locals())
			sj_files = glob.glob('arion/geo_illumina/s05-primates.dir/{organism}/STAR/pass1/*/*-SJ.out.tab'.format(**locals()))
			outfile = 'arion/geo_illumina/s05-primates.dir/{organism}/STAR/pass2/{sample_name}/{sample_name}-Aligned.sortedByCoord.out.bam'.format(**locals())
			yield [[fastq_pair, star_index, sj_files], outfile]

@follows(buildStarIndex, getPrimateStarJunctions)

@files(primateStarJobs)

def runPrimateStar(infiles, outfile):

	# Split infiles
	fastq_files, star_index, sj_files = infiles

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
	fastq_str = ' '.join(fastq_files)
	sj_files_str = ' '.join(sj_files)

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_str} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 5000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="02:00", GB=10, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_outfile=True, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# find arion/geo_illumina/s05-primates.dir/*/STAR/pass2 -name "*.log" | jsc

#############################################
########## 6. RSEM expression
#############################################

# @follows(runPrimateStar)

@transform('arion/geo_illumina/s05-primates.dir/*/STAR/pass2/*/*-Aligned.toTranscriptome.out.bam',
		   regex(r'(.*)/STAR/.*/(.*)-Aligned.toTranscriptome.out.bam'),
		   add_inputs(r'\1/RSEM/index/*_lifted.idx.fa'),
		   r'\1/RSEM/results/\2/\2.isoforms.results')

def runPrimateRsem(infiles, outfile):

	# Variables
	prefix = outfile[:-len('.isoforms.results')]
	reference_name = infiles[1][:-len('.idx.fa')]

	# Paired end
	dataset = outfile.split('/')[-4]
	paired_str = '' if 'macaque_1C_Rep3' in outfile else '--paired-end'
	strandedness_str = 'none'

	# Command
	cmd_str = '''rsem-calculate-expression \
		--alignments \
		--strandedness {strandedness_str} \
		{paired_str} \
		--estimate-rspd \
		--num-threads 200 \
		--no-bam-output \
		{infiles[0]} \
		{reference_name} \
		{prefix} > {prefix}.rsem.log && \
		rsem-plot-model {prefix} {prefix}.quant.pdf '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="02:00", GB=2, n=25, modules=['rsem/1.3.3'], print_outfile=True, stdout=outfile.replace('.isoforms.results', '.log'), stderr=outfile.replace('.isoforms.results', '.err'))

#############################################
########## 7. Prepare metadata
#############################################

# @collate('arion/geo_illumina/s03-alignment.dir/*/RSEM/*/*.isoforms.results',
@collate(runPrimateRsem,
		 regex(r'(.*)/(.*)/RSEM/.*.isoforms.results'),
		 r'\1/\2/RSEM/counts/\2-sample_metadata.tsv')

def preparePrimateSampleMetadata(infiles, outfile):

	# Create directory
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Prepare dict
	sample_dict = [{'files': x, 'names': x.split('/')[-2], 'cell_type': x.split('/')[-2].split('_')[1]} for x in infiles]

	# Convert to dataframe
	sample_dataframe = pd.DataFrame(sample_dict)

	# Write
	sample_dataframe.to_csv(outfile, sep='\t', index=False)

#############################################
########## 8. Aggregate
#############################################

@transform(preparePrimateSampleMetadata,
# @transform('arion/geo_illumina/s05-primates.dir/*/RSEM/counts/*-sample_metadata.tsv',
		   suffix('-sample_metadata.tsv'),
		   add_inputs(human_filtered_gtf),
		   '-counts.rda')

def aggregatePrimateCounts(infiles, outfile):

	# Run
	run_r_job('aggregate_counts', infiles, outfile, conda_env='env', run_locally=False, W='00:10', GB=15, n=1, wait=False, ow=True, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 6. Filter BAM
#############################################

@transform('arion/geo_illumina/s05-primates.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
		   suffix('.bam'),
		   '.filtered.bam')

def filterPrimateBam(infile, outfile):

	# Files
	outname = outfile.replace('.bam', '')

	# Command
	cmd_str = ''' sambamba view --with-header --nthreads 30 --format bam --filter "not unmapped and not duplicate and not secondary_alignment and mapping_quality >= 30" {infile} > {outfile} && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='03:00', n=5, GB=5, modules=['samtools/1.9', 'sambamba/0.5.6'], print_outfile=True)

#############################################
########## 9. Create scaled BigWig
#############################################

# @transform('arion/geo_illumina/s05-primates.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
# 		   regex(r'(.*.dir)/(.*?)/(.*)-Aligned.sortedByCoord.out.bam'),
# 		   add_inputs(r'arion/geo_illumina/s05-primates.dir/\2/RSEM/counts/\2-size_factors.tsv'),
# 		   r'\1/\2/\3-scaled_3.bw')

# def createPrimateScaledBigWig(infiles, outfile):

# 	# Read size factor
# 	normalization_dict = pd.read_table(infiles[1], index_col='sample_name')['size_factor_reciprocal'].to_dict()
# 	size_factor = normalization_dict[outfile.split('/')[-2]]

# 	# Command
# 	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=3 --skipNonCoveredRegions --numberOfProcessors=48 --scaleFactor {size_factor} -b {infiles[0]} -o {outfile}""".format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', W='03:00', n=8, GB=4, print_cmd=False)

#############################################
########## 6. Create BigWig
#############################################

@transform('arion/geo_illumina/s05-primates.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.filtered.bam',
		   suffix('.bam'),
		   '.bw')

def createPrimateScaledBigWig(infile, outfile):

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=3 --normalizeUsing CPM --skipNonCoveredRegions --numberOfProcessors=48 -b {infile} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='03:00', n=8, GB=4, print_outfile=True)

#############################################
########## 7. Merge scaled BigWig
#############################################

# @collate(createScaledBigWig,
# @collate('arion/geo_illumina/s05-primates.dir/*/STAR/pass2/*/*-scaled_3.bw',
@collate('arion/geo_illumina/s05-primates.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.filtered.bw',
		#  regex(r'(.*)/s04-alignment.dir/(.*)/all/STAR/pass2/.*?_(morula)_.*/.*-scaled_3.bw'),
		 regex(r'(.*)/(.*?)/STAR/pass2/(.*?_.*?)_.*/.*.bw'),
		 add_inputs(r'arion/datasets/reference_genomes/\2/*.chrom.sizes'),
		 r'\1/\2/RSEM/counts/filtered_bw/\3.bw')

def mergePrimateScaledBigWig(infiles, outfile):
	
	# Files
	wig_file = outfile.replace('.bw', '.wig')
	bedgraph_file = outfile.replace('.bw', '.bedgraph')
	infiles_str = ' '.join([x[0] for x in infiles])

	# Command
	cmd_str = """ wiggletools mean {infiles_str} > {wig_file} && wigToBigWig {wig_file} {infiles[0][1]} {outfile} && rm {wig_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} | sed '/^KI.*/{{s///;q;}}' > {wig_file} && wigToBigWig {wig_file} {infiles[0][1]} {outfile} && rm {wig_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} > {wig_file} && wiggletools write_bg - {wig_file} | sort -k1,1 -k2,2n > {bedgraph_file} && bedGraphToBigWig {bedgraph_file} {infiles[0][1]} {outfile} && rm {wig_file} {bedgraph_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} | sed '/^KI.*/{{s///;q;}}' > {wig_file} && wiggletools write_bg - {wig_file} | sort -k1,1 -k2,2n > {bedgraph_file} && bedGraphToBigWig {bedgraph_file} {infiles[0][1]} {outfile} && rm {wig_file} {bedgraph_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} > {wig_file} && wiggletools write_bg - {wig_file} | sort -k1,1 -k2,2n > {bedgraph_file} && bedGraphToBigWig {bedgraph_file} {infiles[0][1]} {outfile} && rm {wig_file} {bedgraph_file} """.format(**locals())
	# cmd_str = """ wiggletools mean {infiles_str} > {wig_file} """.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['wiggletools/1.2', 'ucsc-utils/2020-03-17'], W='00:30', n=1, GB=10, print_outfile=False, stdout=outfile.replace('.bw', '.log'), stderr=outfile.replace('.bw', '.err'))

#######################################################
#######################################################
########## S6. ESC datasets
#######################################################
#######################################################

#############################################
########## 1. Create GTF
#############################################

@transform('arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon.gtf',
		   regex(r'.*/(.*).gtf'),
		   add_inputs('arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.gtf'),
		   r'arion/geo_illumina/s06-esc_datasets/indices/gtf/\1.SJ_novel_filtered.gtf')

def createGTF(infiles, outfile):

	# Run
	run_r_job('create_gtf', infiles, outfile, conda_env='env', modules=[], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 2. STAR index
#############################################

@transform(createGTF,
		   regex(r'(.*)/gtf/.*.gtf'),
		   add_inputs('arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'),
		   r'\1/STAR/')

def buildEscStarIndex(infiles, outfile):

	# Command
	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {infiles[1]} --sjdbGTFfile {infiles[0]} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='02:00', GB=5, n=15, ow=True, print_cmd=False, stdout=os.path.join(outfile, 'job.log'), jobname='_'.join(outfile.split('/')[-4:]), wait=False)

#############################################
########## 3. Run bulk ESC STAR
#############################################

@transform('arion/datasets/hesc/rawdata/*.fastq.gz',
		   regex(r'.*/(.*).fastq.gz'),
		   add_inputs('arion/geo_illumina/s06-esc_datasets/indices/STAR'),
		   r'arion/geo_illumina/s06-esc_datasets/bulk/STAR/\1/\1-Aligned.sortedByCoord.out.bam')

def runEscStar(infiles, outfile):

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {infiles[1]} \
		--readFilesIn {infiles[0]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="03:00", GB=10, n=5, modules=['star/2.7.5b', 'samtools/1.11'], print_cmd=True, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.log'))

# ls /hpc/users/torred23/pipelines/projects/early-embryo/arion/geo_illumina/s03-alignment.dir/asami/STAR/pass2/*/*-Aligned.sortedByCoord.out.strandedness.txt | xargs cat

#############################################
########## 4. Strandedness
#############################################

@transform(runEscStar,
		   suffix('.bam'),
		   add_inputs('arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.102_transcript.nochr.bed'),
		   '.strandedness.txt')

def inferStrandedness(infiles, outfile):

	# Command
	cmd_str = ''' infer_experiment.py -i {infiles[0]} -r {infiles[1]} > {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env=False, modules=['python/3.7.3'], W='00:15', GB=10, n=1, print_outfile=False, print_cmd=False)

#############################################
########## 4. Run bulk ESC TEtranscripts
#############################################

@transform('/hpc/users/torred23/pipelines/projects/early-embryo/arion/geo_illumina/s06-esc_datasets/bulk/STAR/*/*-Aligned.sortedByCoord.out.bam',
		   regex(r'(.*)/(.*)-Aligned.sortedByCoord.out.bam'),
		   add_inputs('/hpc/users/torred23/pipelines/projects/early-embryo/geo-illumina/arion/geo_illumina/s06-esc_datasets/indices/gtf/Homo_sapiens.GRCh38.102_talon.SJ_novel_filtered.gtf', '/sc/arion/projects/GuccioneLab/genome-indices/hg38/repeats/GRCh38_Ensembl_rmsk_TE.gtf'),
		   r'\1/\2-TEcount.cntTable')

def runTECount(infiles, outfile):

	# Get prefix
	outdir = os.path.dirname(outfile)
	prefix = os.path.basename(outfile)[:-len('.cntTable')]

	# Command
	cmd_str = ''' cd {outdir} && TEcount -b {infiles[0]} --GTF {infiles[1]} --TE {infiles[2]} --stranded no --sortByPos --project {prefix} '''.format(**locals()).replace('  ', ' ')

	# Run
	run_job(cmd_str, outfile, conda_env=False, modules=['python/3.7.3'], W='10:00', GB=30, n=1, print_outfile=False, print_cmd=False, stdout=outfile.replace('.cntTable', '.log'), stderr=outfile.replace('.cntTable', '.err'))

# find /hpc/users/torred23/pipelines/projects/early-embryo/arion/geo_illumina/s06-esc_datasets -name "*TEcount*" | xargs rm

#############################################
########## 5. Run bulk ESC RSEM
#############################################

#######################################################
#######################################################
########## S7. Single-cell RNA-Seq data
#######################################################
#######################################################

#############################################
########## 1. Filter GTF
#############################################

@transform('arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.102.gtf',
		   regex(r'.*/(.*).gtf'),
		   add_inputs('arion/isoseq/s05-talon.dir/human/gtf/Homo_sapiens.GRCh38.102_talon-SJ_filtered.gtf', 'arion/isoseq/s05-talon.dir/human/gtf/sqanti/Homo_sapiens.GRCh38.102_talon_classification.txt'),
		   r'arion/geo_illumina/s07-scrna.dir/cellranger/reference/\1_talon_SJ_filtered_cellranger.gtf')

def filterGTF(infiles, outfile):

	# Run
	run_r_job('filter_gtf', infiles, outfile, modules=['R/4.0.3'], W='00:15', GB=15, n=1)

#############################################
########## 3. Make cellranger reference
#############################################

@transform('/hpc/users/torred23/pipelines/projects/early-embryo/arion/geo_illumina/s07-scrna.dir/cellranger/reference/Homo_sapiens.GRCh38.102_talon_SJ_filtered_cellranger.gtf',
		   suffix('.gtf'),
		   add_inputs('/hpc/users/torred23/pipelines/projects/early-embryo/geo-illumina/arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'),
		   '_reference')

def makeCellrangerReference(infiles, outfile):

	# Get directory
	outdir = os.path.dirname(outfile)
	basename = os.path.basename(outfile)

	# Command
	cmd_str = ''' cd {outdir} && cellranger mkref --genome={basename} --fasta={infiles[1]} --genes={infiles[0]} --memgb 50 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['cellranger/7.0.0'], W='10:00', GB=15, n=5, print_cmd=True, stdout=os.path.join(outfile, 'job.log'), stderr=os.path.join(outfile, 'job.err'))

#############################################
########## 4. Make STAR index
#############################################

@transform(filterGTF,
		   regex(r'(.*)/cellranger/.*/.*.gtf'),
		   add_inputs('arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa'),
		   r'\1/STAR/index')

def buildStarSoloIndex(infiles, outfile):

	# Command
	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {infiles[1]} --sjdbGTFfile {infiles[0]} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals())

	# Run
	# run_job(cmd_str, outfile, modules=['star/2.7.9a'], W='02:00', GB=5, n=15, ow=True, print_cmd=True, stdout=os.path.join(outfile, 'job.log'), jobname='_'.join(outfile.split('/')[-3:]), wait=False)

#############################################
########## 6. Run STARsolo SmartSeq
#############################################

@transform('arion/datasets/mazid/mazid-smartseq_sample_manifest.tsv',
		   regex(r'.*/(.*)-.*.tsv'),
		   add_inputs(buildStarSoloIndex),
		   r'arion/geo_illumina/s07-scrna.dir/STAR/\1/\1-Aligned.sortedByCoord.out.bam')

def runStarSoloMazid(infiles, outfile):

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {infiles[1]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 64 \
		--soloType SmartSeq \
		--readFilesManifest {infiles[0]} \
		--soloUMIdedup Exact NoDedup \
		--soloStrand Unstranded \
		--soloMultiMappers EM \
		--limitBAMsortRAM 50000000000 \
		--outSAMattributes NH HI nM AS RG GX GN \
		--outSAMtype BAM Unsorted
	'''.format(**locals()) #SortedByCoordinate
		# --limitBAMsortRAM 375913424631 \

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.9a'], W='20:00', GB=10, n=10, print_outfile=False, print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# find /hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/mazid/rawdata/smartseq -name "*CL100183244_L02_1_2*"

#############################################
########## 6. Run STARsolo SmartSeq
#############################################

@transform('arion/datasets/petropoulos/petropoulos-sample_manifest.tsv',
		   regex(r'.*/(.*)-.*.tsv'),
		   add_inputs(buildStarSoloIndex),
		   r'arion/geo_illumina/s07-scrna.dir/STAR/\1/\1-Aligned.sortedByCoord.out.bam')

def runStarSoloPetropoulos(infiles, outfile):

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {infiles[1]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 64 \
		--soloType SmartSeq \
		--readFilesManifest {infiles[0]} \
		--soloUMIdedup Exact NoDedup \
		--soloStrand Unstranded \
		--soloMultiMappers EM \
		--limitBAMsortRAM 50000000000 \
		--outSAMattributes NH HI nM AS RG GX GN \
		--outSAMtype BAM Unsorted
	'''.format(**locals()) #SortedByCoordinate
		# --limitBAMsortRAM 375913424631 \

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.9a'], W='20:00', GB=10, n=10, print_outfile=False, print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# find /hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/mazid/rawdata/smartseq -name "*CL100183244_L02_1_2*"

#############################################
########## 6. Run STARsolo SmartSeq
#############################################

@transform('arion/datasets/yan/yan-sample_manifest.tsv',
		   regex(r'.*/(.*)-.*.tsv'),
		   add_inputs(buildStarSoloIndex),
		   r'arion/geo_illumina/s07-scrna.dir/STAR/\1/\1-Aligned.sortedByCoord.out.bam')

def runStarSoloYan(infiles, outfile):

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {infiles[1]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 64 \
		--soloType SmartSeq \
		--readFilesManifest {infiles[0]} \
		--soloUMIdedup Exact NoDedup \
		--soloStrand Unstranded \
		--soloMultiMappers EM \
		--limitBAMsortRAM 50000000000 \
		--outSAMattributes NH HI nM AS RG GX GN \
		--outSAMtype BAM Unsorted
	'''.format(**locals()) #SortedByCoordinate
		# --limitBAMsortRAM 375913424631 \

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.9a'], W='20:00', GB=10, n=10, print_outfile=False, print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# find /hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/mazid/rawdata/smartseq -name "*CL100183244_L02_1_2*"

#############################################
########## 6. Run STARsolo SmartSeq
#############################################

@transform('arion/datasets/kagawa/kagawa-sample_manifest.tsv',
		   regex(r'.*/(.*)-.*.tsv'),
		   add_inputs(buildStarSoloIndex),
		   r'arion/geo_illumina/s07-scrna.dir/STAR/\1/\1-Aligned.sortedByCoord.out.bam')

def runStarSoloKagawa(infiles, outfile):

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {infiles[1]} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 64 \
		--soloType SmartSeq \
		--readFilesManifest {infiles[0]} \
		--soloUMIdedup Exact NoDedup \
		--soloStrand Unstranded \
		--soloMultiMappers EM \
		--limitBAMsortRAM 50000000000 \
		--outSAMattributes NH HI nM AS RG GX GN \
		--outSAMtype BAM Unsorted
	'''.format(**locals()) #SortedByCoordinate
		# --limitBAMsortRAM 375913424631 \

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.9a'], W='20:00', GB=10, n=10, print_outfile=False, print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# find /hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/mazid/rawdata/smartseq -name "*CL100183244_L02_1_2*"

#############################################
########## 6. Run STARsolo 10X v3
#############################################

@collate(('arion/datasets/*taub*/rawdata/*/*.fastq.gz', 'arion/datasets/liu_iblastoid/merged_lanes/*/*.fastq.gz', 'arion/datasets/yu/merged_lanes/*/*.fastq.gz'),
		 regex(r'.*/datasets/(.*?)/.*?/(.*)/(.*)_S.*_R._.*.fastq.gz'),
		 add_inputs(buildStarSoloIndex, 'arion/datasets/scrna/3M-february-2018.txt'),
		 r'arion/geo_illumina/s07-scrna.dir/STAR/\1/\2/\3-Aligned.sortedByCoord.out.bam')

def runStarSolo10Xv3(infiles, outfile):

	# Variables
	fastq_files = [x[0] for x in infiles]
	fastq_files.sort()
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]

	# Command
	cmd_str = ''' STAR \
		--readFilesIn {fastq_files[1]} {fastq_files[0]} \
		--readFilesCommand zcat \
		--genomeDir {infiles[0][1]} \
		--outFileNamePrefix {prefix} \
		--runThreadN 30 \
		--soloCellFilter EmptyDrops_CR \
		--soloMultiMappers EM \
		--soloType CB_UMI_Simple \
		--soloCBwhitelist {infiles[0][2]} \
		--soloCBlen 16 \
		--soloUMIstart 17 \
		--soloUMIlen 12 \
		--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
		--outSAMtype BAM SortedByCoordinate
	'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['star/2.7.9a'], W='03:00', GB=10, n=5, print_outfile=False, print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# find /hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/mazid/rawdata/smartseq -name "*CL100183244_L02_1_2*"

#############################################
########## 7. Read STARsolo results
#############################################

# arion/geo_illumina/s07-scrna.dir/STAR/*/*/*-Solo.out/Gene/filtered

@transform('arion/geo_illumina/s07-scrna.dir/STAR/*/*-Solo.out/Gene/filtered',
		   regex(r'(.*)/STAR/.*/(.*)-Solo.out/Gene/filtered'),
		   r'\1/seurat/\2/qc/\2-seurat.rda')

# @transform('arion/geo_illumina/s07-scrna.dir/STAR/*/*/*-Solo.out/Gene/filtered',
# 		   regex(r'(.*)/STAR/(.*?)/.*/(.*)-Solo.out/Gene/filtered'),
# 		   r'\1/seurat/\2/qc/\3-seurat.rda')

def createSeuratObject(infile, outfile):

	# Run
	# if not 'old' in outfile:
	run_r_job('create_seurat_object', infile, outfile, modules=['R/4.0.3'], W='00:30', GB=30, n=1, print_outfile=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 8. Filter Seurat data
#############################################

@transform(createSeuratObject,
		#    regex(r'(.*)/qc/(.*).rda'),
		   regex(r'(.*kagawa.*)/qc/(.*).rda'),
		#    r'\1/sctransform/\2_filtered_sctransform.rda')
		   r'\1/normalized/\2_filtered_normalized.rda')

def filterSeuratObject(infile, outfile):

	# Run
	run_r_job('filter_seurat_object', infile, outfile, modules=['geos/3.10.3', 'R/4.0.3'], W='00:30', GB=50, n=1, print_outfile=False, print_cmd=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 8. Plot UMAP
#############################################

def umapJobs():
	for infile in ['arion/geo_illumina/s07-scrna.dir/seurat/kagawa/normalized/kagawa-seurat_filtered_normalized.rda']:
		for n_neighbors in [10, 20, 30, 40, 50, 60]:
			for min_dist in [.1, .2, .3, .4, .5, .6, .7, .8, .9, 1]:
				outfile = infile.replace('normalized', 'umap').replace('.rda', '_umap_{n_neighbors}_{min_dist}.tsv'.format(**locals()))
				yield [infile, outfile, n_neighbors, min_dist]

@files(umapJobs)

def runUmap(infile, outfile, n_neighbors, min_dist):

	# Run
	run_r_job('run_umap', infile, outfile, additional_params=[n_neighbors, min_dist], modules=['R/4.0.3'], W='00:05', GB=10, n=1, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 1. STAR junctions
#############################################

# @follows(linkFASTQ)

@collate('arion/datasets/kagawa/rawdata/smartseq_merged/*.fastq.gz',
		 regex(r'.*/(.*)/rawdata/.*/(.*)_R..fastq.gz'),
		 add_inputs(human_star_index),
		 r'arion/geo_illumina/s07-scrna.dir/STAR_merged/\1/\2/pass1/\1_\2-SJ.out.tab')

def getStarJunctionsSmartseq(infiles, outfile):

	# Split
	fastq_files = [x[0] for x in infiles]
	star_index = infiles[0][1]

	# FASTQ string
	fastq_str = ' '.join(fastq_files)

	# Prefix
	prefix = outfile[:-len('SJ.out.tab')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_str} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 100 \
		--outSAMtype None'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="06:00", print_cmd=False, GB=5, n=15, modules=['star/2.7.5b'], stdout=outfile.replace('-SJ.out.tab', '_job.log'), q='premium')

# ls arion/geo_illumina/s03-alignment.dir/*/STAR/pass1/*/*_job.log | jsc

#############################################
########## 2. STAR BAM
#############################################

# @follows(getStarJunctionsSmartseq)

@collate('arion/datasets/mazid/rawdata/smartseq_merged/*.fastq.gz',
		 regex(r'.*/(.*)/rawdata/.*/(.*)_R..fastq.gz'),
		 add_inputs(human_star_index, r'arion/geo_illumina/s07-scrna.dir/STAR_merged/\1/\2/pass1/*-SJ.out.tab'),
		 r'arion/geo_illumina/s07-scrna.dir/STAR_merged/\1/\2/pass2/\1_\2-Aligned.sortedByCoord.out.bam')

def runStarSmartseq(infiles, outfile):

	# Split
	fastq_files = [x[0] for x in infiles]
	star_index = infiles[0][1]
	sj_files = infiles[0][2:]

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
	fastq_str = ' '.join(fastq_files)
	sj_files_str = ' '.join(sj_files)

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_str} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 5000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM Unsorted '''.format(**locals())
		# --outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="20:00", GB=30, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_outfile=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# find arion/geo_illumina/s07-scrna.dir/STAR_merged -name "*.log" | js
# find arion/geo_illumina/s07-scrna.dir/STAR_merged -name "*.err" | lr

#############################################
########## 1. STAR junctions
#############################################

# @follows(linkFASTQ)

@transform('arion/datasets/petropoulos/rawdata/smartseq_merged/*.fastq.gz',
		   regex(r'.*/(.*)/rawdata/.*/(.*).fastq.gz'),
		   add_inputs(human_star_index),
		   r'arion/geo_illumina/s07-scrna.dir/STAR_merged/\1/\2/pass1/\1_\2-SJ.out.tab')

def getStarJunctionsSmartseqPetropoulos(infiles, outfile):

	# Split
	fastq_file = infiles[0]
	star_index = infiles[1]

	# Prefix
	prefix = outfile[:-len('SJ.out.tab')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_file} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 100 \
		--outSAMtype None'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="06:00", print_cmd=False, GB=5, n=15, modules=['star/2.7.5b'], stdout=outfile.replace('-SJ.out.tab', '_job.log'), q='premium')

# ls arion/geo_illumina/s03-alignment.dir/*/STAR/pass1/*/*_job.log | jsc

#############################################
########## 2. STAR BAM
#############################################

# @follows(getStarJunctionsSmartseq)

@transform('arion/datasets/petropoulos/rawdata/smartseq_merged/*.fastq.gz',
		   regex(r'.*/(.*)/rawdata/.*/(.*).fastq.gz'),
		   add_inputs(human_star_index, r'arion/geo_illumina/s07-scrna.dir/STAR_merged/\1/\2/pass1/*-SJ.out.tab'),
		   r'arion/geo_illumina/s07-scrna.dir/STAR_merged/\1/\2/pass2/\1_\2-Aligned.sortedByCoord.out.bam')

def runStarSmartseqPetropoulos(infiles, outfile):

	# Split
	fastq_file = infiles[0]
	star_index = infiles[1]
	sj_files = infiles[2:]

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
	sj_files_str = ' '.join(sj_files)

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_file} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 5000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="30:00", GB=15, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))


#############################################
########## 2. STAR BAM
#############################################

def mergeKagawaJobs():
	cluster_dataframe = pd.read_table('arion/geo_illumina/s07-scrna.dir/seurat/kagawa/normalized/kagawa-seurat_filtered_normalized_clusters.tsv')
	cluster_dataframe['R1'] = [glob.glob('arion/datasets/kagawa/rawdata/smartseq/{x}/*R1*.fastq.gz'.format(**locals()))[0] for x in cluster_dataframe['barcode']]
	cluster_dataframe['R2'] = [glob.glob('arion/datasets/kagawa/rawdata/smartseq/{x}/*R2*.fastq.gz'.format(**locals()))[0] for x in cluster_dataframe['barcode']]
	for cluster, read_dataframe in cluster_dataframe.groupby('cluster'):
		for read_mate in ['R1', 'R2']:
			infiles = read_dataframe[read_mate].tolist()
			outfile = 'arion/datasets/kagawa/rawdata/smartseq_merged/kagawa_{cluster}_{read_mate}.fastq.gz'.format(**locals())
			yield [infiles, outfile]

@files(mergeKagawaJobs)

def mergeKagawaFastq(infiles, outfile):

	# Command
	infiles_str = ' '.join(infiles)
	cmd_str = ''' cat {infiles_str} > {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='00:30', GB=15, n=1, print_cmd=False)

#############################################
########## . Individual STAR
#############################################

@collate(('arion/datasets/kagawa/rawdata/smartseq/*/*_R*.f*.gz', 'arion/datasets/mazid/rawdata/smartseq/*/*.f*.gz'),
# @collate(('arion/datasets/kagawa/rawdata/smartseq/SRR14798089/*_R*.f*.gz', 'arion/datasets/mazid/rawdata/smartseq/CNR0462970/*.fq.gz'),
		 regex(r'.*/(.*)/rawdata/.*/(.*)/.*.gz'),
		 add_inputs(human_star_index, r'arion/geo_illumina/s07-scrna.dir/STAR/\1/\1-SJ.out.tab'),
		 r'arion/geo_illumina/s07-scrna.dir/STAR_individual/\1/\2/\2-Aligned.sortedByCoord.out.bam')

def runStarSmartseqIndividual(infiles, outfile):
	
	# Split
	fastq_files = [x[0] for x in infiles]
	star_index = infiles[0][1]
	sj_files = infiles[0][2:]

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
	fastq_str = ' '.join(fastq_files)
	sj_files_str = ' '.join(sj_files)

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_str} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--sjdbFileChrStartEnd {sj_files_str} \
		--limitSjdbInsertNsj 5000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())
		# --outSAMtype BAM Unsorted '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="01:00", GB=10, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# find /hpc/users/torred23/pipelines/projects/early-embryo/arion/geo_illumina/s07-scrna.dir/ -name "*sorted.20h.bam.tmp*" | xargs rm
# ls arion/geo_illumina/s07-scrna.dir/STAR_individual/*/*/*-Log.out | lr > arion/geo_illumina/s07-scrna.dir/STAR_individual/star_logs.txt

#############################################
########## 2. STAR BAM
#############################################

# @follows(getStarJunctionsSmartseq)

@transform('arion/datasets/petropoulos/rawdata/smartseq/*.fastq.gz',
		   regex(r'.*/(.*)/rawdata/.*/(.*).fastq.gz'),
		   add_inputs(human_star_index, r'arion/geo_illumina/s07-scrna.dir/STAR/\1/\1-SJ.out.tab'),
		 r'arion/geo_illumina/s07-scrna.dir/STAR_individual/\1/\2/\2-Aligned.sortedByCoord.out.bam')

def runStarSmartseqPetropoulosIndividual(infiles, outfile):

	# Split
	fastq_file = infiles[0]
	star_index = infiles[1]
	sj_file = infiles[2]

	# Variables
	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]

	# Command
	cmd_str = ''' STAR \
		--genomeDir {star_index} \
		--readFilesIn {fastq_file} \
		--readFilesCommand zcat \
		--outFileNamePrefix {prefix} \
		--runThreadN 32 \
		--sjdbFileChrStartEnd {sj_file} \
		--limitSjdbInsertNsj 5000000 \
		--quantMode TranscriptomeSAM GeneCounts \
		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="01:00", GB=10, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

#############################################
########## 4. RSEM expression
#############################################

@transform('arion/geo_illumina/s07-scrna.dir/STAR_individual/*/*/*-Aligned.toTranscriptome.out.bam',
# @transform(('arion/geo_illumina/s07-scrna.dir/STAR_individual/kagawa/SRR14795911/SRR14795911-Aligned.toTranscriptome.out.bam', 'arion/geo_illumina/s07-scrna.dir/STAR_individual/mazid/CNR0462970/CNR0462970-Aligned.toTranscriptome.out.bam', 'arion/geo_illumina/s07-scrna.dir/STAR_individual/petropoulos/ERR1041403/ERR1041403-Aligned.toTranscriptome.out.bam'),
		   regex(r'(.*)/STAR_individual/(.*)/(.*)/.*-Aligned.toTranscriptome.out.bam'),
		   add_inputs(human_rsem_index),
		   r'\1/RSEM_individual/\2/\3/\3.isoforms.results')

def runRsemSmartseq(infiles, outfile):

	# Variables
	prefix = outfile[:-len('.isoforms.results')]
	reference_name = infiles[1][:-len('.idx.fa')]

	# Paired end
	dataset = outfile.split('/')[-4]
	paired_str = '' if 'petropoulos' in outfile else '--paired-end'

	# Command
	cmd_str = '''rsem-calculate-expression \
		--alignments \
		--strandedness none \
		{paired_str} \
		--estimate-rspd \
		--num-threads 200 \
		--no-bam-output \
		{infiles[0]} \
		{reference_name} \
		{prefix} > {prefix}.rsem.log '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W="03:00", GB=2, n=25, modules=['rsem/1.3.3'], print_cmd=False, stdout=outfile.replace('.isoforms.results', '.log'), stderr=outfile.replace('.isoforms.results', '.err'))

# ls arion/geo_illumina/s07-scrna.dir/RSEM_individual/petropoulos/*/*.isoforms.results | wc -l
# ls arion/geo_illumina/s07-scrna.dir/RSEM_individual/petropoulos/*/*.err | wc -l

# #############################################
# ########## . Merged SJ support
# #############################################

# @transform('arion/geo_illumina/s07-scrna.dir/STAR_merged/mazid/*/pass1/*-SJ.out.tab',
# 		   regex(r'(.*)/pass1/(.*)-SJ.out.tab'),
# 		#    add_inputs(additionalInputs),
# 		   r'\1/SJ_support/\2-SJ_support.txt')

# def getMergedSJSupport(infiles, outfile):
# 	print(infiles, outfile)


# #############################################
# ########## 8. Filter Seurat data
# #############################################

# @collate(createSeuratObject,
# 		 regex(r'(.*)/(.*liu.*)/qc/.*-seurat.rda'),
# 		 r'\1/\2/filtered_himem_1/\2-seurat_filtered.rda')

# def filterSeuratObject(infiles, outfile):

# 	# Run
# 	run_r_job('filter_seurat_object', infiles, outfile, modules=['R/4.0.3'], W='03:00', GB=30, n=5, run_locally=False, print_outfile=False, print_cmd=False, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

# Taubenschmid-Stowers
# Paper - cells with more than 2.500 and less than 10.500 genes, as well as less than 10% mitochondrial reads, and less than 10% largest gene percentage were kept
# Our params - >2,500, <10,000, <10% mito

# Yu Blastoid - min nFeature, max nFeature, max mito
# Paper - Cutoffs of nGene < 1,300 and pctMT > 15 were applied to discard cells
# 36 - 3K, 10K, 10 - tissue: blastoids, internal sample: LW36
# 49 - 2.5K, 8K, 10 - Naïve WIBR3; 5iLA
# 50 - 1.5K, 5K, 10 - Blastoid-nESCs; 5iLA
# 51 - 2K, 6K, 10 - Blastoid-nEND; NACL
# 52 - 1.3K, 5K, 15 (lots of MT reads) - Blastoid -hTSCs; hTSM
# 58 - 800, 3K, 15 (lots of MT) - internal sample: LW58, naive hescs cultrue medium: naive hESCs grown in 5iLA, collection day: D3
# 59 - 1K, 3K, 10 - internal sample: LW59,vnaive hescs cultrue medium: naive hESCs grown in 5iLA, collection day: D6
# 60 - 1K, 4K, 15 - internal sample: LW60, naive hescs cultrue medium: naive hESCs grown in 5iLA, collection day: D9
# 61 - 800, 3K, 15 - internal sample: LW61, naive hescs cultrue medium: naive hESCs grown in PXGL, collection day: D9

# Liu iBlastoid
# Scrublet was used to detect multiplets (v.0.2.1).
# Cells with more than 20% of unique molecular identifiers (UMIs) from mitochondria were filtered out.
# Cells with more than 200 genes detected were included in the analysis. Genes detected in fewer than 30 cells or with fewer than 60 UMIs in total across all single cells were filtered before normalization

# Mazid
# For our data, cells with fewer than 200 genes detected, greater than 12,000 genes detected or more than 10% mitochondrial reads were removed.
# For SMART-seq2 data, cells with fewer than 1,500,000 counts were filtered out from downstream analysis. 
# For DNBelab C4 data of ES cells, cells with fewer than 200 genes detected, greater than 12,000 genes detected or more than 10% mitochondrial reads were filtered out

# Petropoulos
# Cells were quality-filtered based on four criteria, leaving 1,529 cells post-filtering out of 1,919 sequenced cells. 
# First, Spearman correlations, using the RPKM expression levels of all genes, for every possible pair of cells were calculated and a histogram of the
# maximum correlation obtained for each cell, corresponding to the most similar cell, was used to identify 305 outlier cells with a maximum pair-wise correlations below 0.63 (Figure 1A). 
# Second, a histogram of the number of expressed genes per cell was used to identify 330 outlier cells with 
# less than 5000 expressed genes (Figure 1B).
# Third, a histogram of the total transcriptional expression output from the sex chromosomes (RPKM sum)
# was used to identify 33 cells with indeterminable sex, or a called sex that was inconsistent with other cells of that
# embryo (see Figure S1 and section “Calling Embryonic Sex”, below). 
# Fourth, 13 outlier cells were identified using PCA and t-SNE dimensionality reduction. 

#############################################
########## 3. Run CellRanger
#############################################

@transform('/hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/liu_iblastoid/rawdata',
		   regex(r'(.*)/datasets/(.*)/rawdata'),
		   add_inputs('/hpc/users/torred23/pipelines/projects/early-embryo/arion/geo_illumina/s07-scrna.dir/cellranger/reference/Homo_sapiens.GRCh38.102_talon_SJ_filtered_cellranger_reference'),
		   r'\1/geo_illumina/s07-scrna.dir/cellranger/\2/\2-cellranger_results.txt')

def runCellRanger(infiles, outfile):

	# Variables
	outdir = os.path.dirname(outfile)
	basename = os.path.basename(outfile)[:-len('-cellranger_results.txt')]

	# Command
	cmd_str = ''' cd {outdir} && cellranger count \
		--id {basename} \
		--transcriptome {infiles[1]} \
		--fastqs {infiles[0]} \
		--localcores 10
	'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['cellranger/7.0.0'], W='07:00', GB=10, n=10, print_outfile=False, print_cmd=False, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'))

#############################################
########## 3. Run CellRanger
#############################################

@subdivide('/hpc/users/torred23/pipelines/projects/early-embryo/arion/datasets/liu_iblastoid/rawdata',
		   regex(r'(.*)/datasets/(.*)/rawdata'),
		   add_inputs('/hpc/users/torred23/pipelines/projects/early-embryo/arion/geo_illumina/s07-scrna.dir/cellranger/reference/Homo_sapiens.GRCh38.102_talon_SJ_filtered_cellranger_reference'),
		   r'\1/geo_illumina/s07-scrna.dir/cellranger/\2/*/\2-cellranger_results.txt',
		   r'\1/geo_illumina/s07-scrna.dir/cellranger/\2/{sample_name}/\2_{sample_name}-cellranger_results.txt')

def runCellRangerLiu(infiles, outfiles, outfileRoot):

	# Get samples
	sample_list = set([os.path.basename(x).split('_')[0] for x in glob.glob(os.path.join(infiles[0], '*', '*'))])

	# Loop
	for sample_name in sample_list:

		# Variables
		outfile = outfileRoot.format(**locals())
		outdir = os.path.dirname(outfile)
		basename = os.path.basename(outfile)[:-len('-cellranger_results.txt')]+'_'+sample_name

		# Command
		cmd_str = ''' cd {outdir} && cellranger count \
			--id {basename} \
			--transcriptome {infiles[1]} \
			--fastqs {infiles[0]} \
			--sample {sample_name} \
			--localcores 10
		'''.format(**locals())

		# Run
		run_job(cmd_str, outfile, modules=['cellranger/7.0.0'], W='05:00', GB=10, n=10, print_outfile=False, print_cmd=False, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'))

#############################################
########## 3. Sort BAM
#############################################

@transform('arion/geo_illumina/s07-scrna.dir/STAR/*/*-Aligned.out.bam',
		   suffix('.bam'),
		   '.sorted.20h.bam')

def sortBam(infile, outfile):

	# Command
	cmd_str = ''' samtools sort -o {outfile} -@ 10 {infile} && samtools index {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['samtools/1.13'], W='20:00', GB=20, n=5, print_outfile=False, print_cmd=False)

#############################################
########## 3. Sort BAM
#############################################

@transform('arion/geo_illumina/s07-scrna.dir/STAR/*/*/*-Aligned.sortedByCoord.out.bam',
		   suffix('.bam'),
		   '.bam.bai')

def indexBam(infile, outfile):

	# Command
	cmd_str = ''' samtools index -@ 32 {infile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['samtools/1.13'], W='03:00', GB=20, n=5, print_outfile=False, print_cmd=False)

#############################################
########## 3. Create scTE reference
#############################################

# @transform(createGTF,
# 		   regex(r'(.*)/s06.*/(.*).gtf'),
# 		   add_inputs(genome_indices['hg38']['ensembl_TE_gtf']),
# 		   r'\1/s07-scrna.dir/scTE/index/\2_scTE_reference.exclusive.idx')

# def createScteReference(infiles, outfile):

# 	# Get parameters
# 	te_bed = os.path.join(os.path.dirname(outfile), os.path.basename(infiles[1])).replace('.gtf', '.bed')
# 	basename = outfile.replace('.exclusive.idx', '')

# 	# Command gtf2bed < {infiles[1]} > {te_bed} && \
		
# 	cmd_str = ''' gtf2bed < {infiles[1]} > {te_bed} && \
# 		scTE_build \
# 			-te {te_bed} \
# 			-gene {infiles[0]} \
# 			-g hg38 \
# 			-o {basename}
# 	'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', modules=['bedops/2.4.36'], W='03:00', GB=30, n=1, print_outfile=False, print_cmd=False, stdout=outfile.replace('.exclusive.idx', '.log'), stderr=outfile.replace('.exclusive.idx', '.log'))

# # Successfully installed bitstring-3.1.9 scTE-1.1.4

#############################################
########## 4. Run scTE
#############################################

# # @transform(runStarSolo,
# @transform('/hpc/users/torred23/pipelines/projects/early-embryo/arion/geo_illumina/s07-scrna.dir/STAR/*/*-Aligned.sortedByCoord.out.bam',
# 		   regex(r'(.*)/STAR/(.*)/(.*)-Aligned.sortedByCoord.out.bam'),
# 		   add_inputs('/sc/arion/scratch/torred23/scrna/mazid/scte/reference/custom_reference.exclusive.idx'),
# 		   r'\1/scTE/\2/\3-scTE')

# def runscTE(infiles, outfile):

# 	# Get outdir
# 	outdir = os.path.dirname(outfile)
# 	basename = os.path.basename(outfile)

# 	# Command
# 	cmd_str = ''' cd {outdir} && scTE \
# 		-i {infiles[0]} \
# 		-x {infiles[1]} \
# 		-o {basename} \
# 		--hdf5 True \
# 		--genome hg38 \
# 		-p 5 \
# 		-CB True \
# 		-UMI True
# 	'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['python/3.7.3', 'samtools/1.13'], W='10:00', GB=15, n=5, print_outfile=False, print_cmd=False, stdout=outfile+'_job.log', stderr=outfile+'_job.err')

	# --soloType CB_UMI_Complex \
	# --soloUMIlen 10 \
	# --soloCBposition 0_0_0_9 0_16_0_25 \
	# --soloUMIposition 0_31_0_40 \
	# --soloCBwhitelist /sc/arion/scratch/torred23/scrna/mazid/whitelist1.txt /sc/arion/scratch/torred23/scrna/mazid/whitelist2.txt \
	# --soloCBmatchWLtype 1MM \

	# --quantMode TranscriptomeSAM GeneCounts \
# @files(('/sc/arion/scratch/torred23/scrna/mazid/WT-CO-1_L4_137A37.R1.fastq.gz', '/sc/arion/scratch/torred23/scrna/mazid/WT-CO-1_L4_137A37.R2.fastq.gz', 'arion/geo_illumina/s06-esc_datasets/indices/STAR', 'arion/datasets/scrna/3M-february-2018.txt'),
# 		'/sc/arion/scratch/torred23/scrna/mazid/star/test3/test-Aligned.sortedByCoord.out.bam')

# # rm -r /sc/arion/scratch/torred23/scrna/mazid/star/test1

# def runStarSolo(infiles, outfile):

# 	# Variables
# 	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {infiles[2]} \
# 		--readFilesIn {infiles[0]} {infiles[1]} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 32 \
# 		--soloType CB_UMI_Complex \
# 		--soloUMIlen 10 \
#         --soloCBposition 0_0_0_9 0_16_0_25 \
#         --soloUMIposition 0_31_0_40 \
# 		--soloCBwhitelist /sc/arion/scratch/torred23/scrna/mazid/whitelist1.txt /sc/arion/scratch/torred23/scrna/mazid/whitelist2.txt \
# 		--soloCBmatchWLtype 1MM \
# 		--outSAMattributes NH HI nM AS CR UR CB UB GX GN sS sQ sM \
# 		--outSAMtype BAM SortedByCoordinate
# 	'''.format(**locals())
# 		# --soloCBwhitelist {infiles[2]} \
# 		# --soloCBlen 16 \
# 		# --soloUMIstart 17 \
# 		# --soloUMIlen 12 \

# 		# --quantMode TranscriptomeSAM GeneCounts \

# 	# Run
# 	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='03:00', GB=10, n=10, print_outfile=False, print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.log'))

# @files((runStarSolo, '/sc/arion/scratch/torred23/scrna/mazid/scte/reference/custom_reference.exclusive.idx'),
# 		'/sc/arion/scratch/torred23/scrna/mazid/scte/results/scte_results')

# # rm -r /sc/arion/scratch/torred23/scrna/mazid/star/test1

# def runscTE(infiles, outfile):

# 	# Command
# 	cmd_str = ''' scTE \
# 		-i {infiles[0]} \
# 		-x {infiles[1]} \
# 		-o {outfile} \
# 		--hdf5 True \
# 		--genome hg38 \
# 		-p 5 \
# 		-CB True \
# 		-UMI True
# 	'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['python/3.7.3', 'samtools/1.13'], W='10:00', GB=15, n=5, print_outfile=False, print_cmd=False, stdout=outfile+'_job.log', stderr=outfile+'_job.err')

# @transform('arion/datasets/scrna/rawdata/SRR14853531.fastq.gz',
# 		   regex(r'.*/(.*).fastq.gz'),
# 		   add_inputs('arion/geo_illumina/s06-esc_datasets/indices/STAR', 'arion/datasets/scrna/3M-february-2018.txt'),
# 		   r'arion/geo_illumina/s06-esc_datasets/scrna/\1/\1-Aligned.sortedByCoord.out.bam')

# def runStarSolo(infiles, outfile):

# 	# Variables
# 	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {infiles[1]} \
# 		--readFilesIn {infiles[0]} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 100 \
# 		--soloCBwhitelist {infiles[2]} \
# 		--soloCBlen 16 \
# 		--soloUMIstart 17 \
# 		--soloUMIlen 12 \
# 		--soloType CB_UMI_Simple \
# 		--quantMode TranscriptomeSAM GeneCounts \
# 		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 100 '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['star/2.7.5b', 'samtools/1.11'], W='06:00', GB=10, n=10, print_outfile=False, print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.log'))

# https://www.biostars.org/p/462568/
# 10x v3
# Whitelist, 3M-Feb_2018_V3.txt
# CB length, 16
# UMI start, 17
# UMI length, 12

# whitelists are available from: https://github.com/10XGenomics/cellranger/tree/master/lib/python/cellranger/barcodes

# STAR \
# ...
# --soloCBwhitelist [whitelist] \
# --soloCBlen [CB length] \
# --soloUMIstart [UMI start] \
# --soloUMIlen [UMI length] \
# ...

#############################################
########## 7. Run SCTE
#############################################

#############################################
########## 6. Filter BAM
#############################################

@transform(runEscStar,
		   suffix('.bam'),
		   '.filtered.bam')

def filterEscBam(infile, outfile):

	# Files
	outname = outfile.replace('.bam', '')

	# Command  and mapping_quality >= 30
	cmd_str = ''' sambamba view --with-header --nthreads 30 --format bam --filter "not unmapped and not duplicate and not secondary_alignment" {infile} > {outfile} && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='03:00', n=5, GB=5, modules=['samtools/1.9', 'sambamba/0.5.6'], print_outfile=False)

#############################################
########## 8. Merge BAM
#############################################

def mergeJobs():
	sample_dict = pd.read_csv('arion/datasets/hesc/hesc-samples.csv', comment='#').groupby('source_name')['Run'].apply(lambda x: list(x))
	for cell_type, samples in sample_dict.items():
		infiles = ['arion/geo_illumina/s06-esc_datasets/bulk/STAR/{x}/{x}-Aligned.sortedByCoord.out.filtered.bam'.format(**locals()) for x in samples]
		outfile = 'arion/geo_illumina/s06-esc_datasets/bulk/merged_bam/{cell_type}-merged.filtered.bam'.format(**locals())
		yield [infiles, outfile]

@files(mergeJobs)

def mergeEscBam(infiles, outfile):

	# Command
	cmd_str = ''' samtools merge {outfile} {infiles[0]} {infiles[1]} -@ 32 && samtools index {outfile} -@ 32 '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='03:00', n=5, GB=5, modules=['samtools/1.11'], print_cmd=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

#############################################
########## 8. Create bigwig
#############################################

# @transform('arion/geo_illumina/s07-scrna.dir/STAR/liu_iblastoid/RL2157_2020_06_24_20_2175_32F38F1_scRNAseq_S6/RL2157_2020_06_24_20_2175_32F38F1_scRNAseq-Aligned.sortedByCoord.out.bam',
@transform('arion/geo_illumina/s07-scrna.dir/STAR_merged/*/*/pass2/*-Aligned.sortedByCoord.out.bam',
# @transform('arion/geo_illumina/s07-scrna.dir/STAR_merged/mazid/8CLC/pass2/mazid_8CLC-Aligned.sortedByCoord.out.bam',
# @transform('arion/geo_illumina/s07-scrna.dir/STAR/mazid/mazid-Aligned.out.sorted.bam',
# @transform('arion/geo_illumina/s07-scrna.dir/STAR/taubenschmid_stowers/SRR14853531/HNES1_naive_hPSC-Aligned.sortedByCoord.out.bam',
		   suffix('.bam'),
		   '.bw')

def createBigWig(infile, outfile):

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=3 --normalizeUsing CPM --skipNonCoveredRegions --numberOfProcessors=48 -b {infile} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='03:00', n=8, GB=4, print_cmd=False)


# # @transform(runPrimateStar,
# # @transform('arion/geo_illumina/s05-primates.dir/*/STAR/pass2/*/*-Aligned.sortedByCoord.out.bam',
# @transform('arion/geo_illumina/s05-primates.dir/marmoset/STAR/pass2/*/*_1_*-Aligned.sortedByCoord.out.bam',
# 		   suffix('-Aligned.sortedByCoord.out.bam'),
# 		   '.bw')

# def createPrimateBigWig(infile, outfile):

# 	# Command
# 	cmd_str = """bamCoverage --outFileFormat=bigwig --skipNonCoveredRegions --numberOfProcessors=50 --normalizeUsing RPKM -b {infile} -o {outfile}""".format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', W='03:00', n=5, GB=6)

#############################################
########## 7. Aggregate counts
#############################################

# find arion/geo_illumina/s05-primates.dir/*/RSEM/results -name "*.log" | grep -v 'rsem.log' | jsc
# find arion/geo_illumina/s05-primates.dir/*/RSEM/results -name "*.isoforms.results" | wc -l

# @follows(linkFASTQ)

# @collate(macaque_linked_fastq,
# 		 regex(r'(.*)/s01-datasets.dir/.*/fastq/(.*)/.*.fastq.gz'),
# 		 add_inputs(buildMacaqueStarIndex),
# 		 r'\1/s05-macaque.dir/alignment/STAR/pass1/\2/\2-SJ.out.tab')

# def getMacaqueStarJunctions(infiles, outfile):

# 	# Split
# 	fastq_files = [x[0] for x in infiles]
# 	star_index = infiles[0][1]

# 	# FASTQ string
# 	fastq_str = ' '.join(fastq_files)

# 	# Prefix
# 	prefix = outfile[:-len('SJ.out.tab')]

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {star_index} \
# 		--readFilesIn {fastq_str} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 100 \
# 		--outSAMtype None'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="02:00", print_outfile=False, GB=5, n=15, modules=['star/2.7.5b'], stdout=outfile.replace('-SJ.out.tab', '_job.log'))

# # find arion/geo_illumina/s05-macaque.dir -name "*job.log" | jsc
# # find arion/geo_illumina/s05-macaque.dir -name "*final.out" | xargs grep -e 'Uniquely'

# #######################################################
# #######################################################
# ########## S5. Macaque analysis
# #######################################################
# #######################################################

# #############################################
# ########## 1. STAR index
# #############################################

# @files((macaque_primary_assembly, macaque_gtf),
# 	   'arion/geo_illumina/s05-macaque.dir/indices/STAR')

# def buildMacaqueStarIndex(infiles, outfile):

# 	# Command
# 	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {infiles[0]} --sjdbGTFfile {infiles[1]} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='02:00', GB=5, n=15, ow=True, print_cmd=True, stdout=os.path.join(outfile, 'job.log'), jobname='_'.join(outfile.split('/')[-4:]), wait=False)

# #############################################
# ########## 2. RSEM index
# #############################################

# @files((macaque_primary_assembly, macaque_gtf),
# 	   'arion/geo_illumina/s05-macaque.dir/indices/RSEM/Macaca_mulatta.Mmul_10.102.idx.fa')

# def buildMacaqueRsemIndex(infiles, outfile):

# 	# Command
# 	basename = outfile[:-len('.idx.fa')]
# 	cmd_str = ''' rsem-prepare-reference --gtf {infiles[1]} --num-threads 10 {infiles[0]} {basename} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="00:30", GB=5, n=3, modules=['rsem/1.3.3'], print_cmd=True, stdout=basename+'.log', stderr=basename+'.err')

# #############################################
# ########## 3. STAR first pass
# #############################################

# # @follows(linkFASTQ)

# @collate(macaque_linked_fastq,
# 		 regex(r'(.*)/s01-datasets.dir/.*/fastq/(.*)/.*.fastq.gz'),
# 		 add_inputs(buildMacaqueStarIndex),
# 		 r'\1/s05-macaque.dir/alignment/STAR/pass1/\2/\2-SJ.out.tab')

# def getMacaqueStarJunctions(infiles, outfile):

# 	# Split
# 	fastq_files = [x[0] for x in infiles]
# 	star_index = infiles[0][1]

# 	# FASTQ string
# 	fastq_str = ' '.join(fastq_files)

# 	# Prefix
# 	prefix = outfile[:-len('SJ.out.tab')]

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {star_index} \
# 		--readFilesIn {fastq_str} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 100 \
# 		--outSAMtype None'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="02:00", print_outfile=False, GB=5, n=15, modules=['star/2.7.5b'], stdout=outfile.replace('-SJ.out.tab', '_job.log'))

# # find arion/geo_illumina/s05-macaque.dir -name "*job.log" | jsc
# # find arion/geo_illumina/s05-macaque.dir -name "*final.out" | xargs grep -e 'Uniquely'

# #############################################
# ########## 4. STAR second pass
# #############################################

# @collate(macaque_linked_fastq,
# 		 regex(r'(.*)/s01-datasets.dir/.*/fastq/(.*)/.*.fastq.gz'),
# 		 add_inputs(buildMacaqueStarIndex, getMacaqueStarJunctions),
# 		 r'\1/s05-macaque.dir/alignment/STAR/pass2/\2/\2-Aligned.sortedByCoord.out.bam')

# def runMacaqueStar(infiles, outfile):

# 	# Split
# 	fastq_files = [x[0] for x in infiles]
# 	star_index = infiles[0][1]
# 	sj_files = infiles[0][2:]

# 	# Variables
# 	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
# 	fastq_str = ' '.join(fastq_files)
# 	sj_files_str = ' '.join(sj_files)

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {star_index} \
# 		--readFilesIn {fastq_str} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 32 \
# 		--sjdbFileChrStartEnd {sj_files_str} \
# 		--limitSjdbInsertNsj 5000000 \
# 		--quantMode TranscriptomeSAM GeneCounts \
# 		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="02:00", GB=15, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_outfile=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# # find arion/geo_illumina/s05-macaque.dir/alignment/STAR/pass2 -name "*.log" | js

# #############################################
# ########## 5. Create BigWig
# #############################################

# @transform(runMacaqueStar,
# 		   suffix('-Aligned.sortedByCoord.out.bam'),
# 		   '.bw')

# def createMacaqueBigWig(infile, outfile):

# 	# Command
# 	cmd_str = """bamCoverage --outFileFormat=bigwig --skipNonCoveredRegions --numberOfProcessors=50 --normalizeUsing RPKM -b {infile} -o {outfile}""".format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', W='05:00', n=7, GB=6)

# #############################################
# ########## 5. RSEM
# #############################################

# #############################################
# ########## 6. Merge
# #############################################

# #######################################################
# #######################################################
# ########## S6. Marmoset analysis
# #######################################################
# #######################################################

# #############################################
# ########## 1. STAR index
# #############################################

# @files((marmoset_primary_assembly, marmoset_gtf),
# 	   'arion/geo_illumina/s06-marmoset.dir/indices/STAR')

# def buildMarmosetStarIndex(infiles, outfile):

# 	# Command
# 	cmd_str = '''STAR --runMode genomeGenerate --genomeDir {outfile} --genomeFastaFiles {infiles[0]} --sjdbGTFfile {infiles[1]} --runThreadN 100 --outFileNamePrefix {outfile}'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['star/2.7.5b'], W='02:00', GB=5, n=15, ow=True, print_cmd=False, stdout=os.path.join(outfile, 'job.log'), jobname='_'.join(outfile.split('/')[-4:]), wait=False)

#############################################
########## 2. RSEM index
#############################################

# @files((marmoset_primary_assembly, marmoset_gtf),
# 	   'arion/geo_illumina/s06-marmoset.dir/indices/RSEM/Macaca_mulatta.Mmul_10.102.idx.fa')

# def buildMarmosetRsemIndex(infiles, outfile):

# 	# Command
# 	basename = outfile[:-len('.idx.fa')]
# 	cmd_str = ''' rsem-prepare-reference --gtf {infiles[1]} --num-threads 10 {infiles[0]} {basename} '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="00:30", GB=5, n=3, modules=['rsem/1.3.3'], print_cmd=True, stdout=basename+'.log', stderr=basename+'.err')

# #############################################
# ########## 3. STAR first pass
# #############################################

# # @follows(linkFASTQ)

# @collate(macaque_linked_fastq,
# 		 regex(r'(.*)/s01-datasets.dir/.*/fastq/(.*)/.*.fastq.gz'),
# 		 add_inputs(buildMacaqueStarIndex),
# 		 r'\1/s05-macaque.dir/alignment/STAR/pass1/\2/\2-SJ.out.tab')

# def getMacaqueStarJunctions(infiles, outfile):

# 	# Split
# 	fastq_files = [x[0] for x in infiles]
# 	star_index = infiles[0][1]

# 	# FASTQ string
# 	fastq_str = ' '.join(fastq_files)

# 	# Prefix
# 	prefix = outfile[:-len('SJ.out.tab')]

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {star_index} \
# 		--readFilesIn {fastq_str} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 100 \
# 		--outSAMtype None'''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="02:00", print_outfile=False, GB=5, n=15, modules=['star/2.7.5b'], stdout=outfile.replace('-SJ.out.tab', '_job.log'))

# # find arion/geo_illumina/s05-macaque.dir -name "*job.log" | jsc
# # find arion/geo_illumina/s05-macaque.dir -name "*final.out" | xargs grep -e 'Uniquely'

# #############################################
# ########## 4. STAR second pass
# #############################################

# @collate(macaque_linked_fastq,
# 		 regex(r'(.*)/s01-datasets.dir/.*/fastq/(.*)/.*.fastq.gz'),
# 		 add_inputs(buildMacaqueStarIndex, getMacaqueStarJunctions),
# 		 r'\1/s05-macaque.dir/alignment/STAR/pass2/\2/\2-Aligned.sortedByCoord.out.bam')

# def runMacaqueStar(infiles, outfile):

# 	# Split
# 	fastq_files = [x[0] for x in infiles]
# 	star_index = infiles[0][1]
# 	sj_files = infiles[0][2:]

# 	# Variables
# 	prefix = outfile[:-len('Aligned.sortedByCoord.out.bam')]
# 	fastq_str = ' '.join(fastq_files)
# 	sj_files_str = ' '.join(sj_files)

# 	# Command
# 	cmd_str = ''' STAR \
# 		--genomeDir {star_index} \
# 		--readFilesIn {fastq_str} \
# 		--readFilesCommand zcat \
# 		--outFileNamePrefix {prefix} \
# 		--runThreadN 32 \
# 		--sjdbFileChrStartEnd {sj_files_str} \
# 		--limitSjdbInsertNsj 5000000 \
# 		--quantMode TranscriptomeSAM GeneCounts \
# 		--outSAMtype BAM SortedByCoordinate && samtools index {outfile} -@ 32 '''.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, W="02:00", GB=15, n=10, modules=['star/2.7.5b', 'samtools/1.11'], print_outfile=False, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))

# # find arion/geo_illumina/s05-macaque.dir/alignment/STAR/pass2 -name "*.log" | js

# #############################################
# ########## 5. Create BigWig
# #############################################

# @transform(runMacaqueStar,
# 		   suffix('-Aligned.sortedByCoord.out.bam'),
# 		   '.bw')

# def createMacaqueBigWig(infile, outfile):

# 	# Command
# 	cmd_str = """bamCoverage --outFileFormat=bigwig --skipNonCoveredRegions --numberOfProcessors=50 --normalizeUsing RPKM -b {infile} -o {outfile}""".format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', W='05:00', n=7, GB=6)

# #############################################
# ########## 5. RSEM
# #############################################

# #############################################
# ########## 6. Merge
# #############################################

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