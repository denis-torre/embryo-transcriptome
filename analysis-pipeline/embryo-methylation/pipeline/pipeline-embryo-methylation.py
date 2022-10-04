#################################################################
#################################################################
############### Embryo methylation ################
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
r_source = 'pipeline/scripts/embryo-methylation.R'
py_source = 'pipeline/scripts/EmbryoMethylation.py'
P = 'acc_apollo'
q = 'sla'
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
#import EmbryoMethylation as P

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
sample_file = 'arion/datasets/guo/guo-samples.csv'
name_file = 'arion/datasets/guo/guo-sample_names.csv'
all_illumina_fastq = glob.glob('arion/methylation/s01-fastq.dir/human/*/*/*.f*q.gz')

#######################################################
#######################################################
########## S1. FASTQ
#######################################################
#######################################################

#############################################
########## 1. Link
#############################################

def linkJobs():
	sample_dataframe = pd.read_csv(sample_file, comment='#').rename(columns={'GEO_Accession (exp)': 'geo_sample_id'})[['Run', 'geo_sample_id']]
	name_dataframe = pd.read_table(name_file)
	name_dataframe['sample_name'] = [x.replace('RRBS', 'human').replace('-cell', 'C_').replace('ICM', 'ICM_').replace('TE', 'TE_').replace('Zygote', '1C_').replace('MII_Oocyte', 'oocyte_').replace('Morula', 'morula_').replace('Postimplantation_embryo', 'postimplantation_').replace('rep', 'Rep') for x in name_dataframe['sample_name']]
	sample_dict = sample_dataframe.merge(name_dataframe, on='geo_sample_id').set_index('Run')['sample_name'].to_dict()
	fastq_files = glob.glob('arion/datasets/guo/rawdata/*.fastq.gz')
	for fastq_file in fastq_files:
		fastq_name = os.path.basename(fastq_file)
		run = fastq_name[:-len('_1.fastq.gz')]
		sample_name = sample_dict[run]
		infile = os.path.join(os.getcwd(), fastq_file)
		outfile = 'arion/methylation/s01-fastq.dir/human/raw/{sample_name}/{fastq_name}'.format(**locals())
		yield [infile, outfile]

@files(linkJobs)

def linkFASTQ(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)

	# Create
	if not os.path.exists(outfile):
		print(outfile)
		# os.system('mkdir -p {outdir} && ln -s {infile} {outfile}'.format(**locals()))

#############################################
########## 2. Adapter trimming
#############################################

def trimJobs():
	for sample_path in glob.glob('arion/methylation/s01-fastq.dir/*/raw/*'):
		infiles = glob.glob(os.path.join(sample_path, '*'))
		infiles.sort()
		outdir = sample_path.replace('/raw/', '/trimmed/')
		yield [infiles, outdir]

# @follows(linkFASTQ)

@files(trimJobs)

def trimIlluminaAdapters(infiles, outdir):

	# Command
	cmd_str = '''trim_galore --paired --rrbs --cores 6 --output_dir {outdir} {infiles[0]} {infiles[1]}'''.format(**locals())

	# Run
	run_job(cmd_str, outdir, modules=['trim_galore/0.6.6'], W='05:00', GB=6, n=6, print_outfile=False, stdout=os.path.join(outdir, 'job.log'), stderr=os.path.join(outdir, 'job.err'))

#######################################################
#######################################################
########## S2. QC
#######################################################
#######################################################

#############################################
########## 1. FASTQC
#############################################

# @follows(linkFASTQ, trimIlluminaAdapters, trimReads)

@transform(all_illumina_fastq,
		   regex(r'(.*)/s01-fastq.dir/(.*).f.*q.gz'),
		   r'\1/s02-fastqc.dir/\2_fastqc.html')

def runFastQC(infile, outfile):

	# Command
	cmd_str = '''fastqc --outdir=$(dirname {outfile}) {infile}'''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['fastqc/0.11.8'], W='03:00', GB=12, n=1, print_outfile=False)

# find arion/chipseq/s02-fastqc.dir/human -name "*fastqc.html" | wc -l

#############################################
########## 2. MultiQC
#############################################

@follows(runFastQC)

def qcJobs():
	filelist = [
		['arion/methylation/s02-fastqc.dir/human/raw', 'arion/methylation/multiqc/human_fastqc/multiqc_report.html'],
		['arion/methylation/s02-fastqc.dir/human/trimmed', 'arion/methylation/multiqc/human_fastqc_trimmed/multiqc_report.html'],
		['arion/methylation/s02-fastqc.dir/human', 'arion/methylation/multiqc/human_fastqc_all/multiqc_report.html']
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
########## S3. Bismark
#######################################################
#######################################################

#############################################
########## 1. Prepare genome
#############################################

@transform('arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly.fa',
		   regex(r'.*/reference_genomes/(.*)/.*.fa'),
		   r'arion/methylation/s03-bismark.dir/\1/genome/Bisulfite_Genome')

def prepareBismarkGenome(infile, outfile):

	# Directory
	genome_dir = os.path.dirname(outfile)

	# Command
	cmd_str = ''' cp {infile} {genome_dir} && bismark_genome_preparation {genome_dir} --verbose '''.format(**locals()) #--parallel 4

	# Run
	run_job(cmd_str, outfile, print_outfile=True, modules=['bismark/0.22.3', 'bowtie2/2.4.1'], W='06:00', GB=50, n=1, stdout=outfile.replace('.fa', '.log'), stderr=outfile.replace('.fa', '.err'))#, print_outfile=False, print_cmd=False)

#############################################
########## 2. Run Bismark
#############################################

@collate('arion/methylation/s01-fastq.dir/human/trimmed/*/*.fq.gz',
		 regex(r'(.*)/s01-fastq.dir/(.*?)/.*/(.*)/.*.fq.gz'),
		 add_inputs(r'\1/s03-bismark.dir/\2/genome/'),
		 r'\1/s03-bismark.dir/\2/results/\3/\3_pe.bam')

def runBismark(infiles, outfile):

	# Directory
	output_dir = os.path.dirname(outfile)
	temp_dir = os.path.join(output_dir, 'tmp')
	basename = os.path.basename(outfile)[:-len('_pe.bam')]

	# Command
	cmd_str = '''  bismark --genome_folder {infiles[0][1]} -1 {infiles[0][0]} -2 {infiles[1][0]} --output_dir {output_dir} --temp_dir {temp_dir} --basename {basename} '''.format(**locals()) #  --parallel 4

	# Run
	run_job(cmd_str, outfile, print_cmd=False, modules=['bismark/0.22.3', 'bowtie2/2.4.1', 'samtools/1.13'], W='30:00', GB=10, n=5, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'))#, print_outfile=False, print_cmd=False)

#############################################
########## 3. Extract methylation
#############################################

@transform('arion/methylation/s03-bismark.dir/human/results/*/*.bam',
		   regex(r'(.*)/(.*).bam'),
		   add_inputs('arion/methylation/s03-bismark.dir/human/genome/'),
		   r'\1/methylation/\2_splitting_report.txt')

def extractMethylation(infiles, outfile):
	
	# Directory
	output_dir = os.path.dirname(outfile)
	genome_folder = os.path.join(os.getcwd(), infiles[1])

	# Command
	cmd_str = '''  bismark_methylation_extractor --gzip --bedGraph --output {output_dir} {infiles[0]} --cytosine_report --genome_folder {genome_folder} --parallel 10 '''.format(**locals()) # --parallel 4 

	# Run
	run_job(cmd_str, outfile, print_cmd=False, modules=['bismark/0.22.3', 'bowtie2/2.4.1', 'samtools/1.13'], W='30:00', GB=5, n=5, stdout=outfile.replace('.txt', '.log'), stderr=outfile.replace('.txt', '.err'))#, print_outfile=False, print_cmd=False)

#############################################
########## 4. Report links
#############################################

# @follows(runBismark)

@transform(('arion/methylation/s03-bismark.dir/human/results/*/*_report.txt', 'arion/methylation/s03-bismark.dir/human/results/*/methylation/*_report.txt', 'arion/methylation/s03-bismark.dir/human/results/*/*.bam'),
		   regex(r'(.*)/results/.*/(.*)'),
		   r'\1/reports/links/\2')

def linkBismarkReports(infile, outfile):

	# Add path
	infile = os.path.join(os.getcwd(), infile)

	# Create directory
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Create link
	if not os.path.exists(outfile):
		os.system('ln -s {infile} {outfile}'.format(**locals()))

#############################################
########## 5. Individual reports
#############################################

# @follows(runBismark)

@transform('arion/methylation/s03-bismark.dir/human/reports/links',
		   regex(r'(.*)/links'),
		   r'\1/split/')

def createBismarkReport(infile, outfile):

	# Add path
	outfile = os.path.join(os.getcwd(), outfile)

	# Command
	cmd_str = ''' cd {infile} && bismark2report --dir {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, print_cmd=False, ow=True, modules=['bismark/0.22.3'], W='00:15', GB=10, n=1, stdout=os.path.join(outfile, 'job.log'), stderr=os.path.join(outfile, 'job.err'))

#############################################
########## 6. Summary report
#############################################

@collate('arion/methylation/s03-bismark.dir/human/reports/links/*.bam',
		 regex(r'(.*)/links/.*.bam'),
		 r'\1/summary/bismark_summary_report.html')

def createBismarkSummary(infiles, outfile):

	# Parameters
	infiles_str = ' '.join(infiles)
	basename = outfile[:-len('.html')]

	# Command
	cmd_str = ''' bismark2summary {infiles_str} --basename {basename} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, print_cmd=False, ow=False, modules=['bismark/0.22.3'], W='00:15', GB=10, n=1, stdout=os.path.join(outfile, 'job.log'), stderr=os.path.join(outfile, 'job.err'))

#######################################################
#######################################################
########## S4. Methylation estimates
#######################################################
#######################################################

#############################################
########## 1. Tile
#############################################

@transform('arion/methylation/s03-bismark.dir/human/results/*/methylation/*.CpG_report.txt.gz',
		   regex(r'(.*)/s03-bismark.dir/(.*)/results/(.*)/methylation/.*.gz'),
		   r'\1/s04-methylation.dir/\2/tiled_100bp/\3-tiled_methylation.tsv.gz')

def tileMethylation(infile, outfile):

	# Run
	run_r_job('tile_methylation', infile, outfile, modules=['R/4.0.3'], W='03:00', GB=30, n=3, stdout=outfile.replace('.tsv.gz', '.log'), stderr=outfile.replace('.tsv.gz', '.err'))

#############################################
########## 2. Average tiled
#############################################

@collate(tileMethylation,
		 regex(r'(.*)/(.*)/(.*)/(human_.*?)_.*.tsv.gz'),
		 add_inputs(r'arion/datasets/reference_genomes/\2/*.chrom.sizes'),
		 r'\1/\2/\3/bw/\4.bw')

def averageTiledMethylation(infiles, outfile):

	# Split
	chromsize_file = infiles[0][1]
	infiles = [x[0] for x in infiles]

	# Run
	run_r_job('average_tiled_methylation', infiles, outfile, modules=['R/4.0.3'], additional_params=chromsize_file, W='03:00', GB=50, n=1, stdout=outfile.replace('.bw', '.log'), stderr=outfile.replace('.bw', '.err'))

#############################################
########## 3. Average bedgraph
#############################################

# @collate(tileMethylation,
# 		 regex(r'(.*)/(.*)/split/(human_.*?)_.*.tsv.gz'),
# 		 add_inputs(r'arion/datasets/reference_genomes/\2/*.chrom.sizes'),
# 		 r'\1/\2/bw_tiled/\3-tiled_methylation.bw')

# def averageTiledMethylation(infiles, outfile):

# 	# Split
# 	chromsize_file = infiles[0][1]
# 	infiles = [x[0] for x in infiles]

# 	# Run
# 	run_r_job('average_tiled_methylation', infiles, outfile, modules=['R/4.0.3'], additional_params=chromsize_file, W='03:00', GB=50, n=1, stdout=outfile.replace('.bw', '.log'), stderr=outfile.replace('.bw', '.err'))

# rm -r /hpc/users/torred23/pipelines/projects/early-embryo/arion/methylation/s04-methylation.dir/human/bw

# #############################################
# ########## 1. Smoothen
# #############################################

# # @collate('arion/methylation/s03-bismark.dir/human/results/*/methylation/*.bismark.cov.gz', # change to whole genome coverage report
# @collate('arion/methylation/s03-bismark.dir/human/results/*/methylation_v2/*.CpG_report.txt.gz', # change to whole genome coverage report
# 		 regex(r'(.*)/s03-bismark.dir/(.*)/results/(.*?_.*?)_.*.gz'),
# 		 r'\1/s04-methylation.dir/\2/cpg/\3-smoothened_methylation.rda')

# def smoothenMethylation(infiles, outfile):

# 	# Run
# 	run_r_job('smoothen_methylation', infiles, outfile, modules=['R/4.0.3'], W='03:00', GB=50, n=1, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

# #############################################
# ########## 2. TSS profiles
# #############################################

# def tssScoreJobs():
# 	for organism in ['human']:
# 		bed_files = glob.glob('arion/atacseq/s07-tss_scores.dir/{organism}/bed/*.bed'.format(**locals()))
# 		methylation_files = glob.glob('arion/methylation/s04-methylation.dir/{organism}/split/*-smoothened_methylation.rda'.format(**locals()))
# 		for methylation_file in methylation_files:
# 			cell_type = os.path.basename(methylation_file).split('-')[0]
# 			for bed_file in bed_files:
# 				bed_name = os.path.basename(bed_file)[:-len('.bed')]
# 				infiles = [methylation_file, bed_file]
# 				outfile = 'arion/methylation/s05-tss_profiles.dir/{organism}/average_scores/{organism}_{cell_type}-{bed_name}-scores.tsv'.format(**locals())
# 				yield [infiles, outfile]

# # @follows(getTssBed, shuffleTssBed, mergeScaledBigWig)

# @files(tssScoreJobs)

# def getTssScores(infiles, outfile):

# 	# Command
# 	cmd_str = ''' bigWigAverageOverBed {infiles[0]} {infiles[1]} {outfile} '''.format(**locals())
	
# 	# Run
# 	print(infiles, outfile)
# 	# run_job(cmd_str, outfile, modules=['ucsc-utils/2020-03-17'], W='00:10', GB=30, n=1, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

# #############################################
# ########## 1. Smoothen
# #############################################

# @transform('arion/methylation/s03-bismark.dir/human/results/*8C*/methylation/*.bismark.cov.gz', # change to whole genome coverage report
# 		   regex(r'(.*)/s03-bismark.dir/(.*)/results/.*/(.*)_pe.bismark.cov.gz'),
# 		   r'\1/s04-methylation.dir/\2/split/\3-smoothened_methylation.rda')

# def smoothenMethylation(infile, outfile):

# 	# Run
# 	run_r_job('smoothen_methylation', infile, outfile, modules=['R/4.0.3'], W='00:30', GB=5, n=5, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

# #############################################
# ########## 2. Combine
# #############################################

# @collate(smoothenMethylation,
# 		 regex(r'(.*)/(.*)/split/.*.rda'),
# 		 r'\1/\2-smoothened_methylation_merged.rda')

# def mergeSmoothenedMethylation(infiles, outfile):
# 	print(infiles, outfile)

# for coverage plots: read each TSS BED into separate GRanges, use getMeth to get methylation across samples,
# average across replicates, merge and plot

# #############################################
# ########## 1. Average replicates
# #############################################

# @collate('arion/methylation/s03-bismark.dir/human/results/*/methylation/*.bedGraph.gz',
# 		 regex(r'(.*)/s03-bismark.dir/(.*)/results/(human_.*?_.).*/.*.gz'),
# 		 r'\1/s04-methylation.dir/\2/bedgraph_split/\3.bedgraph.gz')

# def mergeReplicates(infiles, outfile):

# 	# Check
# 	if len(infiles) == 1:

# 		# Create directory
# 		outdir = os.path.dirname(outfile)
# 		if not os.path.exists(outdir):
# 			os.makedirs(outdir)

# 		# Copy
# 		os.system('cp {infiles[0]} {outfile}'.format(**locals()))

# 	else:

# 		# Run
# 		run_r_job('get_average_methylation', infiles, outfile, conda_env='env', W='00:30', GB=10, n=1, stdout=outfile.replace('.bedgraph.gz', '.log'), stderr=outfile.replace('.bedgraph.gz', '.err'))#, run_locally=False, print_outfile=False, print_cmd=False)

# # find arion/methylation/s04-methylation.dir/human/bedgraph_split -name "*.log" | jsc
# # find arion/methylation/s04-methylation.dir/human/bedgraph_split -name "*.err" | xargs wc -l

# #############################################
# ########## 2. Average samples
# #############################################

# @collate(mergeReplicates,
# 		 regex(r'(.*)/bedgraph_split/(.*)_..bedgraph.gz'),
# 		 r'\1/bedgraph/\2.bedgraph.gz')

# def mergeSamples(infiles, outfile):

# 	# Run
# 	run_r_job('get_average_methylation', infiles, outfile, conda_env='env', W='00:30', GB=10, n=1, stdout=outfile.replace('.bedgraph.gz', '.log'), stderr=outfile.replace('.bedgraph.gz', '.err'))#, run_locally=False, print_outfile=False, print_cmd=False)

# # find arion/methylation/s04-methylation.dir/human/bedgraph -name "*.log" | jsc
# # find arion/methylation/s04-methylation.dir/human/bedgraph -name "*.err" | xargs wc -l

# #############################################
# ########## 3. BigWig
# #############################################

# @transform(mergeSamples,
# 		   regex(r'(.*)/bedgraph/(.*).bedgraph.gz'),
# 		   r'\1/bw/\2.bw')

# def convertBigWig(infile, outfile):

# 	# Temp
# 	tempfile = outfile.replace('.bw', '.sorted.tmp.bedgraph')

# 	# Command
# 	cmd_str = """ zcat {infile} | sort -k1,1 -k2,2n > {tempfile} && bedGraphToBigWig {tempfile} arion/datasets/reference_genomes/human/Homo_sapiens.GRCh38.dna_sm.primary_assembly_renamed.nochr.chromsizes {outfile} && rm {tempfile} """.format(**locals())

# 	# Run
# 	run_job(cmd_str, outfile, modules=['ucsc-utils/2020-03-17'], W='00:30', n=1, GB=10, print_cmd=False, stdout=outfile.replace('.bw', '.log'), stderr=outfile.replace('.bw', '.err'))

# find arion/methylation/s04-methylation.dir/human/bw -name "*.log" | jsc
# find arion/methylation/s04-methylation.dir/human/bw -name "*.err" | xargs wc -l

#######################################################
#######################################################
########## S5. TSS Average
#######################################################
#######################################################

#############################################
########## 1. Get TSS content
#############################################

@transform('arion/atacseq/s07-tss_scores.dir/*/bed*/*bp.bed',
		   regex(r'.*.dir/(.*)/.*/(.*).bed'),
		   add_inputs(r'arion/datasets/reference_genomes/\1/*.dna_sm.primary_assembly.fa'),
		   r'arion/methylation/s05-tss_scores.dir/\1/nucleotide_content/\2-nucleotide_content.tsv')

def getTssContent(infiles, outfile):

	# Command
	cmd_str = ''' bedtools nuc -fi {infiles[1]} -bed {infiles[0]} > {outfile} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env=False, modules=['bedtools/2.29.0'], W='00:15', GB=10, n=1, print_outfile=False, print_cmd=False)

#############################################
########## 3. Get average scores
#############################################

def tssScoreJobs():
	for organism in ['human']:
		bed_files = glob.glob('arion/atacseq/s07-tss_scores.dir/{organism}/bed*/*bp.bed'.format(**locals()))
		cytosine_files = glob.glob('arion/methylation/s03-bismark.dir/{organism}/results/*/methylation/*.CpG_report.txt.gz'.format(**locals()))
		for cytosine_file in cytosine_files:
			sample_name = os.path.basename(cytosine_file).rsplit('_', 2)[0]
			for bed_file in bed_files:
				bed_name = os.path.basename(bed_file)[:-len('.bed')]
				infiles = [cytosine_file, bed_file]
				outfile = 'arion/methylation/s05-tss_scores.dir/{organism}/split/{sample_name}-{bed_name}-scores.tsv'.format(**locals())
				yield [infiles, outfile]

@files(tssScoreJobs)

def getTssMethylation(infiles, outfile):

	# Run
	run_r_job('get_tss_methylation', infiles, outfile, modules=['R/4.0.3'], W='00:30', GB=20, n=1, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#######################################################
#######################################################
########## S5. Plots
#######################################################
#######################################################

#############################################
########## 1. TSS Matrix
#############################################

@follows(createBigWig, splitGTF)

@collate('arion/methylation/s04-methylation.dir/human/tiled_100bp/bw/*.bw',
		 regex(r'(.*)/s04-methylation.dir/(.*)/tiled_100bp/bw/(.*).bw'),
		 add_inputs(r'arion/atacseq/s07-tss_scores.dir/\2/bed*/*_TSS.bed'),
		 r'\1/s06-tss_plots.dir/tss_heatmaps_split/\3/\3-methylation_matrix.gz')

def computeTssMatrixAverage(infiles, outfile):

	# Get order
	order_dict = {
		'bigwig': {'oocyte':0, '1C': 1, '2C': 2, '4C': 3, '8C': 4, 'morula': 5, 'ICM': 6, 'TE': 7, 'postimplantation': 8},
		'bed': {'Known_TSS': 1, 'Novel_TSS': 2, 'Antisense_TSS': 3, 'Intergenic_TSS': 4, 'Shuffled_TSS': 5}
	}

	# Split
	bigwigs = [x[0] for x in infiles] # if 'TBE' not in x[0] and 'control' not in x[0]
	beds = infiles[0][1:]

	# Get bigwig order
	bigwig_str = ' '.join(pd.DataFrame({
		'bigwig': bigwigs,
		'order': [order_dict['bigwig'][os.path.basename(x).split('_')[1].split('.')[0]] for x in bigwigs]
	}).sort_values('order')['bigwig'])

	# Get GTF order
	bed_str = ' '.join(pd.DataFrame({
		'bed': beds,
		'order': [order_dict['bed'][os.path.basename(x).split('.')[0]] for x in beds]
	}).sort_values('order')['bed'])

	# Command
	cmd_str = ''' computeMatrix reference-point -S {bigwig_str} \
					-R {bed_str} \
					--referencePoint center \
					--beforeRegionStartLength 5000 \
					--afterRegionStartLength 5000 \
					--numberOfProcessors 48 \
					-o {outfile}
	'''.format(**locals()) #--skipZeros \

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='06:00', n=6, GB=4, print_cmd=False, ow=False, wait=False)

	# Write samples
	bigwig_names = [os.path.basename(x)[:-len('.bw')].replace('human_', '') for x in bigwig_str.split(' ')]
	jsonfile = outfile.replace('.gz', '.json')
	if not os.path.exists(jsonfile):
		with open(jsonfile, 'w') as openfile:
			openfile.write(json.dumps(bigwig_names))

#############################################
########## 2. TSS Plot
#############################################

@transform('arion/methylation/s06-tss_plots.dir/tss_heatmaps_split/*/*-methylation_matrix.gz',
		   regex(r'(.*).gz'),
		   add_inputs(r'\1.json'),
		   r'\1_profile.png')

def plotTssProfile(infiles, outfile):

	# Read JSON
	with open(infiles[1]) as openfile:
		samples_label = '" "'.join(json.load(openfile))
	data_outfile = outfile.replace('.png', '.txt')

	# Command
	cmd_str = ''' plotProfile -m {infiles[0]} \
					--refPointLabel TSS \
					--samplesLabel "{samples_label}" \
					--plotHeight 6 \
					--plotWidth 6 \
					--colors "#6BAED6" "#EE6A50" "#66C2A4" "#E9967A" "#999999" \
					--legendLocation "lower-right" \
					-out {outfile} \
					--outFileNameData {data_outfile}
	'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='00:30', n=2, GB=5, print_cmd=False, ow=True, run_locally=False)

#############################################
########## 2. Isoform Matrix
#############################################

# @follows(createBigWig, splitGTF)

@collate('arion/methylation/s04-methylation.dir/human/bw_tiled/*.bw',
		 regex(r'(.*)/s04-methylation.dir/(.*)/bw_tiled/.*.bw'),
		 add_inputs(r'arion/atacseq/summary_plots.dir/isoform_heatmaps/\2/gtf/*.gtf'),
		 r'\1/s06-tss_plots.dir/isoform_heatmaps/\2/methylation-matrix.gz')

def computeIsoformMatrixAverage(infiles, outfile):

	# Get order
	order_dict = {
		'bigwig': {'oocyte': 0, '1C': 1, '2C': 2, '4C': 3, '8C': 4, 'morula': 5, 'ICM': 6, 'TE': 7, 'postimplantation': 8},
		'gtf': {'Known': 1, 'NIC': 2, 'NNC': 3, 'Antisense': 4, 'Intergenic': 5}
	}

	# Split
	bigwigs = [x[0] for x in infiles] # if 'TBE' not in x[0] and 'control' not in x[0]
	gtfs = infiles[0][1:]

	# Get bigwig order
	bigwig_str = ' '.join(pd.DataFrame({
		'bigwig': bigwigs,
		'order': [order_dict['bigwig'][os.path.basename(x).replace('-', '_').split('_')[1]] for x in bigwigs]
	}).sort_values('order')['bigwig'])

	# Get GTF order
	gtf_str = ' '.join(pd.DataFrame({
		'gtf': gtfs,
		'order': [order_dict['gtf'][os.path.basename(x).split('.')[0]] for x in gtfs]
	}).sort_values('order')['gtf'])

	# Command
	cmd_str = ''' computeMatrix scale-regions -S {bigwig_str} \
					-R {gtf_str} \
					--beforeRegionStartLength 5000 \
					--regionBodyLength 10000 \
					--afterRegionStartLength 5000 \
					--numberOfProcessors 48 \
					-o {outfile}
	'''.format(**locals()) #--skipZeros 

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='06:00', n=6, GB=4, print_cmd=False, ow=False, wait=False)

	# Write samples
	bigwig_names = [os.path.basename(x)[:-len('.bw')].replace('human_', '') for x in bigwig_str.split(' ')]
	jsonfile = outfile.replace('.gz', '.json')
	if not os.path.exists(jsonfile):
		with open(jsonfile, 'w') as openfile:
			openfile.write(json.dumps(bigwig_names))

#############################################
########## 4. TSS Plot
#############################################

@transform('arion/methylation/s06-tss_plots.dir/isoform_heatmaps_5000/human/methylation-matrix.gz',
		   regex(r'(.*).gz'),
		   add_inputs(r'\1.json'),
		   r'\1_profile.png')

def plotIsoformProfile(infiles, outfile):

	# Read JSON
	with open(infiles[1]) as openfile:
		samples_label = '" "'.join(json.load(openfile))
	data_outfile = outfile.replace('.png', '.txt')

	# Command
	cmd_str = ''' plotProfile -m {infiles[0]} \
					--samplesLabel "{samples_label}" \
					--plotHeight 6 \
					--plotWidth 6 \
					--colors "#6BAED6" "#EE6A50" "#66C2A4" "#E9967A" "#999999" \
					--legendLocation "lower-right" \
					-out {outfile} \
					--outFileNameData {data_outfile}
	'''.format(**locals())
	# --yMax 3 2.7 9 8 16 \

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='00:30', n=2, GB=5, print_cmd=False, ow=True, run_locally=False)

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