#################################################################
#################################################################
############### Embryo ChIP-Seq ################
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
r_source = 'pipeline/scripts/embryo-chipseq.R'
py_source = 'pipeline/scripts/EmbryoChipseq.py'
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
#import EmbryoChipseq as P

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
# All FASTQ
all_illumina_fastq = glob.glob('arion/chipseq/s01-fastq.dir/*/*/*/*.f*q.gz')

# References
reference_dict = {
	'human': {
		'filtered_gtf': 'arion/illumina/s04-alignment.dir/human/all/gtf/Homo_sapiens.GRCh38.102_talon-all-SJ_filtered.gtf',
		'talon_abundance': 'arion/isoseq/s05-talon.dir/human/Homo_sapiens.GRCh38.102_talon_abundance_filtered.tsv'		
	}
}

#######################################################
#######################################################
########## S1. FASTQ
#######################################################
#######################################################

#############################################
########## 1. Link
#############################################

def linkJobs():
	fastq_files = glob.glob('arion/datasets/xia/rawdata/*/*.fastq.gz')
	for fastq_file in fastq_files:
		sample_name = fastq_file.split('/')[-2]
		fastq_name = os.path.basename(fastq_file)
		read_mate = '_'+fastq_name.split('_')[-1][0] if '_' in fastq_name else ''
		infile = os.path.join(os.getcwd(), fastq_file)
		outfile = 'arion/chipseq/s01-fastq.dir/human/raw/{sample_name}/{sample_name}{read_mate}.fastq.gz'.format(**locals())
		yield [infile, outfile]

@files(linkJobs)

def linkFASTQ(infile, outfile):

	# Outdir
	outdir = os.path.dirname(outfile)

	# Create
	if not os.path.exists(outfile):
		os.system('mkdir -p {outdir} && ln -s {infile} {outfile}'.format(**locals()))

#############################################
########## 2. Adapter trimming
#############################################

def trimJobs():
	for sample_path in glob.glob('arion/chipseq/s01-fastq.dir/*/raw/*'):
		infiles = glob.glob(os.path.join(sample_path, '*'))
		infiles.sort()
		outdir = sample_path.replace('/raw/', '/trimmed/')
		yield [infiles, outdir]

# @follows(linkFASTQ)

@files(trimJobs)

def trimIlluminaAdapters(infiles, outdir):

	# Command
	if len(infiles) == 1:
		cmd_str = '''trim_galore --cores 6 --output_dir {outdir} {infiles[0]}'''.format(**locals())
	elif len(infiles) == 2:
		cmd_str = '''trim_galore --paired --cores 6 --output_dir {outdir} {infiles[0]} {infiles[1]}'''.format(**locals())

	# Run
	run_job(cmd_str, outdir, modules=['trim_galore/0.6.6'], W='10:00', GB=6, n=6, print_outfile=True, stdout=os.path.join(outdir, 'job.log'), stderr=os.path.join(outdir, 'job.err'))

#############################################
########## 3. Read trimming
#############################################

def trim50Jobs():
	for sample_path in glob.glob('arion/chipseq/s01-fastq.dir/*/raw/*'):
		infiles = glob.glob(os.path.join(sample_path, '*'))
		infiles.sort()
		outdir = sample_path.replace('/raw/', '/trimmed_50bp/')
		yield [infiles, outdir]

# @follows(linkFASTQ)

@files(trim50Jobs)

def trimReads(infiles, outdir):

	# Command
	cmd_str = '''trim_galore --cores 6 --hardtrim5 50 --output_dir {outdir} {infiles[0]}'''.format(**locals())

	# Run
	run_job(cmd_str, outdir, modules=['trim_galore/0.6.6'], W='06:00', GB=6, n=6, print_cmd=False, stdout=os.path.join(outdir, 'job.log'), stderr=os.path.join(outdir, 'job.err'))

#######################################################
#######################################################
########## S2. QC
#######################################################
#######################################################

#############################################
########## 1. FASTQC
#############################################

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

# @follows(runFastQC)

def qcJobs():
	filelist = [
		['arion/chipseq/s02-fastqc.dir/human/raw', 'arion/chipseq/multiqc/human_fastqc/multiqc_report.html'],
		['arion/chipseq/s02-fastqc.dir/human/trimmed', 'arion/chipseq/multiqc/human_fastqc_trimmed/multiqc_report.html'],
		['arion/chipseq/s02-fastqc.dir/human/trimmed_50bp', 'arion/chipseq/multiqc/human_fastqc_trimmed_50bp/multiqc_report.html'],
		['arion/chipseq/s03-alignment.dir/human', 'arion/chipseq/multiqc/human_alignment_trimmed_50bp/multiqc_report.html']
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
########## S3. Alignment
#######################################################
#######################################################

#############################################
########## 1. Bowtie
#############################################

@follows(trimIlluminaAdapters)

@collate('arion/chipseq/s01-fastq.dir/human/trimmed_50bp/*/*.fq.gz',
 		 regex(r'(.*)/s01-fastq.dir/(.*)/trimmed_50bp/(.*)/.*.fq.gz'),
		 add_inputs(r'arion/chipseq/s03-alignment.dir/\2/bowtie2/index/*primary_assembly.1.bt2'),
		 r'\1/s03-alignment.dir/\2/bowtie2/results/\3/\3.bam')

def runBowtie(infiles, outfile):

	# Infiles
	fastq = [x[0] for x in infiles]
	bowtie_index = infiles[0][1].replace('.1.bt2', '')
	outname = outfile.replace('.bam', '')

	# FASTQ string
	if len(fastq) == 1:
		fastq_str = '-U {fastq[0]} -N 1 -L 25'.format(**locals())
	elif len(fastq) == 2:
		fastq_str = '-1 {fastq[0]} -2 {fastq[1]} -X 1000 --no-mixed --no-discordant'.format(**locals())

	# Command
	cmd_str = '''bowtie2 -x {bowtie_index} {fastq_str} -N 1 -q -p 40 | samtools view -bS --threads 40 | samtools sort --threads 40 -o {outfile} && \
		samtools index -b {outfile} && samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | cut -f5 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \\t]*//' > {outname}.mapq'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='06:00', n=10, GB=5, modules=['bowtie2/2.4.1', 'samtools/1.9'], stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'), print_outfile=False, ow=False)

#############################################
########## 2. Filter
#############################################

@transform(runBowtie,
		   suffix('.bam'),
		   '_filtered.bam')

def filterBam(infile, outfile):

	# Files
	outname = outfile.replace('.bam', '')

	# Command
	cmd_str = ''' sambamba view --with-header --nthreads 30 --format bam --filter "ref_id <= 24 and ref_id != 22 and not unmapped and not duplicate and mapping_quality >= 30" {infile} > {outfile} && \
		samtools index -b {outfile} && samtools flagstat {outfile} > {outname}.flagstat && samtools idxstats {outfile} > {outname}.idxstats && samtools stats {outfile} > {outname}.stats && samtools view {outfile} | awk '$9>0' | cut -f9 | sort | uniq -c | sort -b -k2,2n | sed -e 's/^[ \\t]*//' > {outname}.fragment_lengths.txt'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, W='03:00', n=5, GB=5, modules=['samtools/1.9', 'sambamba/0.5.6'])#, stdout=outfile.replace('.bam', '.log'), stderr=outfile.replace('.bam', '.err'), print_outfile=False)

#############################################
########## 3. Create BigWig
#############################################

# @transform(filterBam,
@transform('arion/chipseq/s03-alignment.dir/human/bowtie2/results/*/*_filtered.bam',
		   suffix('.bam'),
		   add_inputs('/sc/arion/projects/GuccioneLab/genome-indices/hg38/blacklists/hg38-blacklist.v2.bed'),
		   '.bw')

def createBigWig(infiles, outfile):

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=10 --skipNonCoveredRegions --numberOfProcessors=48 --normalizeUsing RPKM --blackListFileName {infiles[1]} -b {infiles[0]} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='06:00', n=8, GB=4, print_outfile=False)

#######################################################
#######################################################
########## S4. Peaks
#######################################################
#######################################################

#############################################
########## 1. MACS2
#############################################

@collate('arion/chipseq/s03-alignment.dir/human/bowtie2/results/*/*_filtered.bam',
		  regex(r'(.*)/s03-alignment.dir/(.*)/bowtie2/results/(.*)_Rep.*/.*.bam'),
		  r'\1/s04-peaks.dir/\2/macs2/\3/\3_peaks.xls')

def runMacs2(infiles, outfile):

	# Base settings
	bam_str = ' '.join(infiles)
	basename = outfile.replace('_peaks.xls', '')
	organism = outfile.split('/')[-4].replace('human', 'hs').replace('mouse', 'mm')

	# Specific settings
	additional_parameters = '--broad --broad-cutoff 0.05'# if 'H3K27me3' in outfile or 'H3K27ac' in outfile else ''
	
	# --nomodel \
	# --nolambda \
	# Command
	cmd_str = ''' macs2 callpeak \
		-t {bam_str} \
		-f BAM \
		-q 0.05 \
		{additional_parameters} \
		-g {organism} \
		-n {basename} '''.format(**locals())

	# Run
	run_job(cmd_str, outfile, modules=['macs/2.1.0'], W='06:00', n=1, GB=30, print_cmd=False, stdout=outfile.replace('.xls', '.log'), stderr=outfile.replace('.xls', '.err'))

#############################################
########## 2. Rename
#############################################

@transform(runMacs2,
		   suffix('.xls'),
		   '.chr.xls')

def renamePeaks(infile, outfile):

	# Read peaks
	peak_dataframe = pd.read_table(infile, comment='#')

	# Add chromosome
	peak_dataframe['chr'] = ['chr{x}'.format(**locals()) for x in peak_dataframe['chr']]
	peak_dataframe['name'] = [os.path.basename(x) for x in peak_dataframe['name']]

	# Write
	peak_dataframe.to_csv(outfile, sep='\t', index=False)

#######################################################
#######################################################
########## S5. Peak Counts
#######################################################
#######################################################

#############################################
########## 1. Sample dataframe
#############################################

@collate('arion/chipseq/s03-alignment.dir/human/bowtie2/results/*/*_filtered.bam',
		 regex(r'(.*)/s03-alignment.dir/(.*)/bowtie2/results/.*_(H3.*?)_Rep.*.bam'),
		 r'\1/s05-counts.dir/\2/\3/\2-\3-chipseq_samples.csv')

def getSampleMetadata(infiles, outfile):

	# Sample dataframe
	sample_dataframe = pd.DataFrame([{
			'SampleID': x.split('/')[-2],
			'bamReads': x,
			'PeakCaller': 'macs',
			'PeakFormat': 'macs'
		} for x in infiles])

	# Add information
	sample_dataframe['Tissue'] = [x.split('_')[1].replace('control', 'control_8C').replace('TBE', 'TBE_8C') for x in sample_dataframe['SampleID']]
	sample_dataframe['Condition'] = [x.split('_H3')[0].split('_')[-1] for x in sample_dataframe['SampleID']]
	sample_dataframe['Factor'] = [x.split('_Rep')[0].split('_')[-1] for x in sample_dataframe['SampleID']]
	sample_dataframe['Peaks'] = ['arion/chipseq/s04-peaks.dir/human/macs2/human_{Tissue}_{Condition}_{Factor}/human_{Tissue}_{Condition}_{Factor}_peaks.chr.xls'.format(**rowData) for index, rowData in sample_dataframe.iterrows()]

	# Outdir
	outdir = os.path.dirname(outfile)
	if not os.path.exists(outdir):
		os.makedirs(outdir)

	# Write
	# pd.set_option('max.colwidth', -1)
	# print(sample_dataframe)
	sample_dataframe.to_csv(outfile, index=False)

#############################################
########## 2. Get counts
#############################################

@transform(getSampleMetadata,
		   suffix('_samples.csv'),
		   '_peak_counts_broad.rda')

def getPeakCounts(infile, outfile):

	# Run
	run_r_job('get_peak_counts', infile, outfile, conda_env='env', W='06:00', GB=100, n=1, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'))

#############################################
########## 3. Get size factors
#############################################

@transform(getPeakCounts,
		   suffix('_peak_counts.rda'),
		   '_size_factors.tsv')

def getSizeFactors(infile, outfile):

	# Run
	run_r_job('get_size_factors', infile, outfile, conda_env='env', W='00:30', GB=50, n=1, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 4. Create scaled BigWig
#############################################

@transform(filterBam,
		   regex(r'(.*)/(s03-alignment.dir)/(.*?)/(.*_)(H3.*?)(_.*).bam'),
		   add_inputs('/sc/arion/projects/GuccioneLab/genome-indices/hg38/blacklists/hg38-blacklist.v2.bed', r'\1/s05-counts.dir/\3/\5/\3-\5-chipseq_size_factors.tsv'),
		   r'\1/\2/\3/\4\5\6_scaled.bw')

def createScaledBigWig(infiles, outfile):

	# Read size factor
	normalization_dict = pd.read_table(infiles[2], index_col='sample_name')['size_factor_reciprocal'].to_dict()
	size_factor = normalization_dict[outfile.split('/')[-2]]

	# Command
	cmd_str = """bamCoverage --outFileFormat=bigwig --binSize=10 --skipNonCoveredRegions --numberOfProcessors=48 --scaleFactor {size_factor} --blackListFileName {infiles[1]} -b {infiles[0]} -o {outfile}""".format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='06:00', n=8, GB=4, print_cmd=False)

#############################################
########## 5. Merge scaled BigWig
#############################################

@collate(createScaledBigWig,
		 regex(r'(.*)/s03-alignment.dir/(.*)/bowtie2/results/(.*)_Rep.*/.*.bw'),
		 r'\1/s05-counts.dir/\2/merged_bw/\3.bw')

def mergeScaledBigWig(infiles, outfile):
	
	# Multiple files
	if len(infiles) > 1:

		# Files
		wig_file = outfile.replace('.bw', '.wig')
		bedgraph_file = outfile.replace('.bw', '.bedgraph')
		infiles_str = ' '.join(infiles)

		# Command
		cmd_str = """ wiggletools mean {infiles_str} > {wig_file} && wiggletools write_bg - {wig_file} > {bedgraph_file} && bedGraphToBigWig {bedgraph_file} arion/chipseq/s05-counts.dir/human/hg38.chrom.sizes {outfile} && rm {wig_file} {bedgraph_file} """.format(**locals())

		# Run
		run_job(cmd_str, outfile, modules=['wiggletools/1.2', 'ucsc-utils/2020-03-17'], W='00:30', n=1, GB=10, print_cmd=False, stdout=outfile.replace('.bw', '.log'), stderr=outfile.replace('.bw', '.err'))

	# Single file
	else:
		if not os.path.exists(outfile):
			os.system('cp {infiles[0]} {outfile}'.format(**locals()))

#############################################
########## 6. Consensus peaks
#############################################

# @transform(getPeakCounts,
@transform('arion/chipseq/s05-counts.dir/human/*/human-*-chipseq_peak_counts*.rda',
		   regex(r'(.*)_peak_counts(.*).rda'),
		   r'\1_consensus_peaks\2.bed')

def getConsensusPeaks(infile, outfile):

	# Run
	run_r_job('get_consensus_peaks', infile, outfile, conda_env='env', W='00:30', GB=10, n=1, stdout=outfile.replace('.bed', '.log'), stderr=outfile.replace('.bed', '.err'))

# find /hpc/users/torred23/pipelines/projects/early-embryo/arion/chipseq/s05-counts.dir/human -name "*consensus*" | xargs rm

#############################################
########## 7. Differential peaks
#############################################

# @transform(getPeakCounts,
@transform(('arion/chipseq/s05-counts.dir/human/H3K4me3/human-H3K4me3-chipseq_peak_counts.rda', 'arion/chipseq/s05-counts.dir/human/H3K27ac/human-H3K27ac-chipseq_peak_counts.rda', 'arion/chipseq/s05-counts.dir/human/H3K27ac/human-H3K27ac-chipseq_peak_counts_broad.rda'),
		   regex(r'(.*)_peak_counts(.*)'),
		   r'\1_differential_peaks\2')

def getDifferentialPeaks(infile, outfile):

	# Run
	run_r_job('get_differential_peaks', infile, outfile, modules=['R/4.0.3'], W='01:00', GB=50, n=1, stdout=outfile.replace('.rda', '.log'), stderr=outfile.replace('.rda', '.err'), print_outfile=False)

#######################################################
#######################################################
########## S6. TSS coverage
#######################################################
#######################################################

#############################################
########## 1. Get TSS BED
#############################################

# def transcriptJobs():
# 	for organism, reference_info in reference_dict.items():
# 		infile = reference_info['filtered_gtf']
# 		outfile = 'arion/chipseq/s06-tss_coverage.dir/{organism}/{organism}-tss.bed'.format(**locals())
# 		yield [infile, outfile]

# @files(transcriptJobs)

# def getTssBed(infile, outfile):

# 	# Run
# 	run_r_job('get_tss_bed', infile, outfile, run_locally=False, conda_env='env')#, modules=[], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 2. Intersect peaks
#############################################

# @follows(runGenrich, getTssRange)

@transform('arion/chipseq/s04-peaks.dir/human/macs2/*/*.broadPeak',
		   regex(r'(.*)/s04-peaks.dir/(.*)/macs2/(.*)/.*.broadPeak'),
		   add_inputs(r'arion/atacseq/s06-tss_coverage.dir/\2/\2-tss_500bp.bed'),
		   r'\1/s06-tss_coverage.dir/\2/intersect/\3-tss_peaks_intersect.bed')

def intersectTssPeaks(infiles, outfile):

	# Command
	cmd_str = ''' bedtools intersect -wa -a {infiles[1]} -b {infiles[0]} > {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['bedtools/2.29.2'], W='00:30', n=1, GB=25, print_cmd=False)#, stdout=outfile.replace('.narrowPeak', '.log'), stderr=outfile.replace('.narrowPeak', '.err'))

#######################################################
#######################################################
########## S6. TSS scores
#######################################################
#######################################################

#############################################
########## 1.1 Split TSS by isoform class
#############################################

# def tssJobs():
# 	for organism, reference_info in reference_dict.items():
# 		infiles = list(reference_info.values())
# 		outfile = 'arion/chipseq/s07-tss_scores.dir/{organism}/bed/Known_TSS.bed'.format(**locals())
# 		yield [infiles, outfile]

# @files(tssJobs)

# def splitTssTypes(infiles, outfile):

# 	# Run
# 	run_r_job('split_tss_types', infiles, outfile, run_locally=False, conda_env='env')#, modules=[], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 2. Shuffle TSS BED
#############################################

# @follows(splitTssTypes)

# @transform('arion/chipseq/s07-tss_scores.dir/*/bed/*.bed',
# 		   regex(r'(.*)/(.*)/bed/(.*).bed'),
# 		   add_inputs(r'arion/datasets/reference_genomes/\2/*.nochr.chromsizes', r'arion/datasets/reference_genomes/\2/*_transcript.bed'),
# 		   r'\1/\2/bed/\3_shuffled.bed')

# def shuffleTssBed(infiles, outfile):

# 	# Command
# 	cmd_str = ''' bedtools shuffle -excl {infiles[2]} -i {infiles[0]} -g {infiles[1]} > {outfile} '''.format(**locals()) # -chrom
	
# 	# Run
# 	run_job(cmd_str, outfile, modules=['bedtools/2.29.2'], W='00:05', GB=5, n=1, print_cmd=False, ow=False)

#############################################
########## 3. TSS overlap
#############################################

def tssScoreJobs():
	for organism in ['human']:
		bed_files = glob.glob('arion/atacseq/s07-tss_scores.dir/{organism}/bed*/*bp.bed'.format(**locals()))
		bigwig_files = glob.glob('arion/chipseq/s05-counts.dir/{organism}/merged_bw/*.bw'.format(**locals()))
		for bigwig_file in bigwig_files:
			bigwig_name = os.path.basename(bigwig_file)[:-len('.bw')]
			for bed_file in bed_files:
				bed_name = os.path.basename(bed_file)[:-len('.bed')]
				infiles = [bigwig_file, bed_file]
				outfile = 'arion/chipseq/s07-tss_scores.dir/{organism}/average_scores/{bigwig_name}-{bed_name}-scores.tsv'.format(**locals())
				yield [infiles, outfile]

# @follows(getTssBed, shuffleTssBed, mergeScaledBigWig)

@files(tssScoreJobs)

def getTssScores(infiles, outfile):

	# Command
	cmd_str = ''' bigWigAverageOverBed {infiles[0]} {infiles[1]} {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['ucsc-utils/2020-03-17'], W='00:10', GB=30, n=1, print_cmd=False, stdout=outfile.replace('.tsv', '.log'), stderr=outfile.replace('.tsv', '.err'))

#############################################
########## 4. Intersect peaks
#############################################

# @follows(runGenrich, getTssRange)
def tssTypeJobs():
	for dataset in ['H3K4me3', 'H3K27ac']:
		peak_file = 'arion/chipseq/s05-counts.dir/human/H3K4me3/human-H3K4me3-chipseq_consensus_peaks_broad.bed' if dataset == 'H3K4me3' else 'arion/chipseq/s05-counts.dir/human/H3K27ac/human-H3K27ac-chipseq_consensus_peaks.bed'
		for tss_type in ['Known', 'Novel', 'Antisense', 'Intergenic']:
			infiles = [peak_file, 'arion/atacseq/s07-tss_scores.dir/human/bed/{tss_type}_TSS.bed'.format(**locals())]
			outfile = 'arion/chipseq/s07-tss_scores.dir/human/tss_consensus_peak_intersection/{dataset}-{tss_type}_consensus_TSS_peak_intersection.bed'.format(**locals())
			yield [infiles, outfile]

@files(tssTypeJobs)

def intersectTssTypePeaks(infiles, outfile):

	# Command
	cmd_str = ''' bedtools intersect -wa -a {infiles[0]} -b {infiles[1]} > {outfile} '''.format(**locals())
	
	# Run
	run_job(cmd_str, outfile, modules=['bedtools/2.29.2'], W='00:30', n=1, GB=25, print_cmd=False)#, stdout=outfile.replace('.narrowPeak', '.log'), stderr=outfile.replace('.narrowPeak', '.err'))

#######################################################
#######################################################
########## Plots
#######################################################
#######################################################

#############################################
#############################################
########## 2. TSS heatmap
#############################################
#############################################

#############################################
########## 1.1 Split TSS by isoform class
#############################################

# def tssJobs():
# 	for organism, reference_info in reference_dict.items():
# 		infiles = list(reference_info.values())
# 		outfile = 'arion/chipseq/summary_plots.dir/tss_heatmaps/{organism}/bed/Known_TSS.bed'.format(**locals())
# 		yield [infiles, outfile]

# @files(tssJobs)

# def splitTSS(infiles, outfile):

# 	# Run
# 	run_r_job('split_tss', infiles, outfile, run_locally=False, conda_env='env')#, modules=[], W='00:15', GB=10, n=1, run_locally=False, print_outfile=False, print_cmd=False)

#############################################
########## 1.2 Get background
#############################################

# @collate('arion/chipseq/s07-tss_scores.dir/human/bed/*_shuffled.bed',
# 		 regex(r'(.*)/s07-tss_scores.dir/(.*)/bed/.*.bed'),
# 		 r'\1/summary_plots.dir/tss_heatmaps/\2/bed/Shuffled_TSS.bed')

# def mergeShuffledTSS(infiles, outfile):

# 	# Get infiles
# 	infiles_str = ' '.join(infiles)

# 	# Merge
# 	os.system('cat {infiles_str} > {outfile}'.format(**locals()))

#############################################
########## 1.3 Matrix
#############################################

# @follows(createBigWig, splitGTF)

@collate('arion/chipseq/s05-counts.dir/human/merged_bw/*.bw',
		 regex(r'(.*)/s05-counts.dir/(.*)/merged_bw/.*_(.*).bw'),
		 add_inputs(r'arion/atacseq/s07-tss_scores.dir/\2/bed*/*_TSS.bed'),
		 r'\1/summary_plots.dir/tss_heatmaps/\2/\3/\2-\3-matrix.gz')

def computeTssMatrixAverage(infiles, outfile):

	# Get order
	order_dict = {
		'bigwig': {'1C': 1, '2C': 2, '4C': 3, '8C': 4, 'morula': 5, 'ICM': 6, 'TE': 7},
		'bed': {'Known_TSS': 1, 'Novel_TSS': 2, 'Antisense_TSS': 3, 'Intergenic_TSS': 4, 'Shuffled_TSS': 5}
	}

	# Split
	bigwigs = [x[0] for x in infiles if 'TBE' not in x[0] and 'control' not in x[0]]
	beds = infiles[0][1:]

	# Get bigwig order
	bigwig_str = ' '.join(pd.DataFrame({
		'bigwig': bigwigs,
		'order': [order_dict['bigwig'][os.path.basename(x).split('_')[1]] for x in bigwigs]
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
					--beforeRegionStartLength 500 \
					--afterRegionStartLength 500 \
					--numberOfProcessors 48 \
					--skipZeros -o {outfile}
	'''.format(**locals())

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='10:00', n=6, GB=4, print_cmd=False, ow=False, wait=False)

	# Write samples
	bigwig_names = [os.path.basename(x)[:-len('.bw')].replace('human_', '').replace('_', ' ').split(' H3')[0] for x in bigwig_str.split(' ')]
	jsonfile = outfile.replace('.gz', '.json')
	if not os.path.exists(jsonfile):
		with open(jsonfile, 'w') as openfile:
			openfile.write(json.dumps(bigwig_names))

#############################################
########## 1.4 Plot
#############################################

@transform(computeTssMatrixAverage,
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
	# --yMax 3 2.7 9 8 16 \

	# Run
	run_job(cmd_str, outfile, conda_env='env', W='00:30', n=2, GB=5, print_cmd=False, ow=True, run_locally=False)

#############################################
########## 1.4 Plot
#############################################

# # 'arion/chipseq/summary_plots.dir/tss_heatmaps_all_reps_scaled/human/human-matrix.gz', 
# @transform(('arion/chipseq/summary_plots.dir/tss_heatmaps/human/human-matrix.gz'),
# # @transform(computeTssMatrix,
# 		   regex(r'(.*).gz'),
# 		   add_inputs(r'\1.json'),
# 		   r'\1_v4.png')

# def plotTssHeatmap(infiles, outfile):

# 	# Read JSON
# 	with open(infiles[1]) as openfile:
# 		samples_label = '" "'.join(json.load(openfile))

# 	# Command
# 					# --samplesLabel "{samples_label}" \
# 	cmd_str = ''' plotHeatmap -m {infiles[0]} \
# 					--heatmapWidth 5 \
# 					--heatmapHeight 10 \
# 					--colorMap Reds \
# 					--missingDataColor 1 \
# 					--refPointLabel TSS \
# 					--legendLocation none \
# 					--zMax 1 \
# 					-out {outfile}
# 	'''.format(**locals())
# 					# --zMax 0.5 1.5 2.5 3 3 2 2.5 \
# 					# --yMin 0 \
# 					# --zMax 60 30 80 120 150 50 60 \

# 	# Run
# 	run_job(cmd_str, outfile, conda_env='env', W='00:30', n=2, GB=5, print_cmd=False, ow=True, run_locally=True)

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