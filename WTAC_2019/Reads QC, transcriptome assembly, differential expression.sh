# Inside the "final_indexes_mm10" folder we can find: there are the fasta sequences of the mouse genome along the 
# ERCC_sequence folder: The fasta file and the gtf file of the ERCC sequences.
# SIRV_sequence: The fasta file of the isoform sequences that correspond to each one of the spliced SIRV sequences.We are not using at all this file. it is only provided for reference.
# SIRV_original_index: The fasta file of the unspliced SIRV sequences. The folder also has the gtf file with the coordinates of each one of the isoforms projected on the unspliced sequence of each one SIRV
# mm10_genome_ens75_and_gtf: The folder has the fasta file of the mouse genome (Mus_musculus.GRCm38.75.dna.toplevel.fa) and the correspoding gene models with (mm10_ens75_SIRV_original_ERCC.gtf) or without (Mus_musculus.GRCm38.75.gtf) the SIRV and ERCC gene models. Obviously, the SIRV and ERCC gene models are not part of mouse genome that will be used latter
# mm10_transcriptome_ens75: The folder has the cDNA of the mm10 transccriptome 

#In the software below we will use the mouse genome and transcriptome to align the nanopore reads 
#Here we will use the mouse genome mm10 with the sequences from the SIRV and ERCC and the corresponding gene models file.
#The non-repeat masked mouse genome mm10 (Ensembl version 96) is downloaded from the Enslembl website (ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/dna/Mus_musculus.GRCm38.dna.toplevel.fa.gz). We add on the fasta file of the mouse genome the fasta sequences of the SIRV and ERCC files downloaded from the manufacturer,s website (https://www.lexogen.com/sirvs/download/).
#The corresponding gene models for the mouse genome are downloaded from the Enslembl website (ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/Mus_musculus.GRCm38.96.gtf.gz). The corresponding isoform models for the SIRV sequences are downloaded from the manufacturer's website (https://www.lexogen.com/sirvs/download/). 
# The mouse transcriptome data are downloaded from ensembl(ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/Mus_musculus.GRCm38.cdna.all.fa.gz)  

#To avoid extensive computation time we will align not in all chromosomes but only on the chr19 genome (mm10_chr19_ens75_SIRV_original_ERCC and 
#mm10_chr19_ens75_SIRV_original_ERCC.gtf)
# we will also use the transcriptome

#From the "final_indexes_mm10" folder on the desktop .
#I used the mm10_genome_ens75_SIRV_original_ERCC as the "genome_of_interest.fa"
#I also used the  mm10_ens75_SIRV_original_ERCC.gtf as the "genome_of_interest.gtf"
#and the corresponding minimap2 and gmap index


#Initially we will start with the sequencing run statistics derived from the summary file or from the fastq file
```
cd summary_file_folder #go inside the "summary file" folder and process the "sequencing_summary.txt" file
```
#The output will be inside the "nanoplot_of_summary_file" folder  
```
NanoPlot -o nanoplot_of_summary_file --readtype 1D --summary sequencing_summary.txt --loglength
```
```
cd ..
```
or
cd sequenced_files #go inside the "sequenced_files" folder and process the "ass_reads.fastq" file
#The output will be inside the "nanoplot_of_the_original_fastq_file" folder 
NanoPlot -o nanoplot_of_the_original_fastq_file -t 2 --fastq pass_reads.fastq  --maxlength 40000 --plots hex dot
cd ..

#OPTIONAL# The pass reads used have an average quality of 7 and up. The reads in the fastq file can be filtered more based on average quality

cd sequenced_files
sequence_filter.py -q 8 pass_reads.fastq above_q8_pass_reads.fastq
cd ..


=====================================================================================================================================


#Here we will use Primer-chop to identify full length reads based on the presence or absence of the adaptor sequences in the beggining or the end of the reads. Afterwards the adaptors will be trimmed from the sequences.
#Primer-chop will create the following files:
#good-fwd.fa: sequences from forward (sense) reads with both primers in the expected positions. The TSO oligo is found at the 5' end. The polyT oligo is found at the 3' end.The adaptors that were found have been trimmed from the sequences.
#good-rev.fa: sequences from reverse (antisense) reads with both primers in the expected positions.The polyT oligo is found at the 5' end. The TSO oligo is found at the 3' end.The adaptors that were found have been trimmed from the sequences. The reverse complement of the original reads is in the files. 
#no_head-fwd.fa: sequences from forward (sense) reads with good tail primer but no detectable head primer. Only the polyT oligo is found at the 3' end.The adaptors that were found have been trimmed from the sequences.
#no_head-rev.fa: sequences from reverse (antisense) reads with good tail primer but no detectable head primer. Only the TSO oligo is found at the 3' end.The adaptors that were found have been trimmed from the sequences.The reverse complement of the original reads is in the files. 
#no_tail-fwd.fa: sequences from forward (sense) reads with good head primer but no detectable tail primer. Only the TSO oligo is found at the 5' end.The adaptors that were found have been trimmed from the sequences.
#no_tail-rev.fa: sequences from reverse (antisense) reads with good head primer but no detectable tail primer. Only the polyT oligo is found at the 3' end.The adaptors that were found have been trimmed from the sequences.The reverse complement of the original reads is in the files. 
#log.txt: IDs of reads with no detected primers.

cd sequenced_files
./../primer-chop/primer-chop -P2 ../primers.fa above_q8_pass_reads.fastq primer_chop_output

#Here we will try to find how many reads have been reported in each one of the produced files. For this we count the number of lines in each of the produced files. As As these are fasta files the number of reads in each file is the number of lines devided by 2.
 
wc -l primer_chop_output/good-fwd.fa
wc -l primer_chop_output/good-rev.fa
wc -l primer_chop_output/no_head-fwd.fa
wc -l primer_chop_output/no_head-rev.fa
wc -l primer_chop_output/no_tail-fwd.fa
wc -l primer_chop_output/no_tail-rev.fa

#take the files that contain  either the  sence or antisense reads that have adaptors present at both end of the molecules  
cat primer_chop_output/good-fwd.fa primer_chop_output/good-rev.fa >adapter_removed_full_length_above_q8_pass_reads.fa

#Because the "pipeline-pinfish-analysis" that we will use in the next step takes only fastq files we will make the fasta file into a dummy fastq file with the reformat.sh script from  the  BBMap suite
#sudo apt install openjdk-11-jre-headless
PATH="$PATH:/home/wtac/Desktop/new_approach/BBMap_38.50b/bbmap/"
reformat.sh in=adapter_removed_full_length_above_q8_pass_reads.fa out=adapter_removed_full_length_above_q8_pass_reads.fastq qfake=12



cd ..

#=====================================================================================================================================

# afterwards we align on the genome in order to identify new spliced isoforms  
# for this we will use the "pipeline-pinfish-analysis" from nanoporetech 
# for processing speed we will only align on 1 chromosome and on all the SIRVs
# 


#add minimap2 in path!!!
PATH="$PATH:/home/wtac/Desktop/new_approach/s1_WT/minimap2-2.17_x64-linux"

cd pipeline_pinfish_analysis_folder #GO INTO THIS FOLDER

#-------

#All the parameters that are needed for the run and that can be changed are in in the config.yml file:
## General pipeline parameters:
# Name of the pipeline:
# pipeline: "pipeline-pinfish-analysis" --> name of the pipeline
# workdir_top: "Workspaces"   --> this is the folder inside the pinfish_folder where the output will be produced
# repo: "https://github.com/nanoporetech/pipeline-pinfish-analysis.git" --> leave it as is
# genome_fasta: "/home/wtac/Desktop/new_approach/s1_WT/chromosomes/chr1.fa" --> Path of the genome in fasta format to be used 
# reads_fastq: "/home/wtac/Desktop/new_approach/s1_WT/quality_8_filtered_pass_reads.fastq" --> Path to the filtered fastq reads to be used

#-------

# minimap_index_opts: "-k14" --> Extra option passed to minimap2 when generating index. Leave it as default
# minimap2_opts: "-uf"--> Extra options passed to minimap2 . The -uf forces minimap2 to consider the forward transcript strand only which makes sense as our nanopore reads they all correspond to the mRNA sequence.
#minimap2_opts: "--splice-flank=no" --> Enable this if you are interested in accuaret mapping of the SIRV data otherwise ignore it (see the minimap2 webpage for explanation).
# minimum_mapping_quality: 10  -->Minmum mapping quality. Leave it as default. Mapping quality is defined as MapQ = -10*log10(P) where P is the probability that this mapping is NOT the correct one for example if it maps to 2 locations then P = 0.5, 3 locations implies P = 2/3, 4 locations => P = 3/4 ... n locations P=n/(n+1) . So MapQ=3 if it maps to 2 locations in the target, MapQ=2 if it maps to 3 locations in the target, MapQ=1  if it maps to 4-9 locations, MapQ=0 if it maps to 10 or more locations and MapQ=255 for unique mapping

#spliced_bam2gff_opts: "-s" -->Options passed to spliced_bam2gff.  A tool for converting sorted BAM files containing spliced alignments (generated by minimap2 or GMAP) into GFF2 format. Each read will be represented as a distinct transcript.

#-------

#The tool takes a sorted GFF2 file as input and clusters together reads having similar exon/intron structure and creates a rough consensus of the clusters by taking the median of exon boundaries from all transcripts in the cluster.
# minimum_cluster_size: 10 --> number of different transcripts in a cluster.Transcript clusters having size less than this parameter are discarded. This parameter has the largest effect on the sensitivity and specificity of transcript reconstruction. Larger values usually lead to higher specificity at the expense of lowering sensitivity.

# Minimum isoform percent
# minimum_isoform_percent: 1.0  --> Fraction of the total reads

# The following parameters are the maximum distance (in bp) tolerated at the start of the first exon and the end of last exon (terminal_exon_boundary_tolerance) , while the Internal exon boundary tolerance is the tolerance for all other exon boundaries.
# exon_boundary_tolerance: 10 -->Internal exon boundary tolerance
# terminal_exon_boundary_tolerance: 30 -->Terminal exon boundary tolerance

#-----

#After defining the clusters then for each cluster  we create an error corrected read by mapping all reads on the read with the median length (using minimap2) and polishing it using racon. The polished reads are mapped back to the genome using minimap2.
#   Extra options for polished reads:
#   minimap2_opts_polished: "-uf"  -->see previously
#   minimap2_opts_polished: "--splice-flank=no"  --> see previously
#   spliced_bam2gff_opts_pol: "-s" --> see previously

#-----

#At the end the GFFs generated by the polished clusters are filtered out from transcripts which are likely to be based on RNA degradation products from the 5' end. The tool clusters the input transcripts into "loci" by the 3' ends and discards transcripts which have a compatible transcripts in the loci with more exons.

# Options passed to collapse_partials when collapsing fragmentation artifacts in clustered and polished transcripts:
#The "collapse_internal_tol" parameter is the exon boundary difference tolerated at internal splice sites, while "collapse_five_tol" and "collapse_three_tol" are the tolerance values at the 5' and 3' end respectively. 
#   collapse_internal_tol: 5 --> Internal exon boundary tolerance for collapsing
#   collapse_five_tol: 5000 --> Five prime boundary tolerance for collapsing
#   collapse_three_tol: 30 ---> Three prime boundary tolerance for collapsing
#   threads: 50 --> Number of threads 

#-----

snakemake --use-conda -j 1 all

#The above command creates the following files:
#Into the "pipeline_pinfish_analysis_folder" inside the "Workspaces/pipeline-pinfish-analysis" folders there are the following folders : 
#in the alignments folder:
#WE WILL USE ON IGV THE FILE BELOW:
#reads_aln_sorted.bam - the input reads aligned to the input genome by minimap2 in BAM format.
#polished_reads_aln_sorted.bam - The spliced alignment of the polished transcripts to the input genome.

# in the results folder:

#original alignment coordinates of reads
#WE WILL USE ON IGV THE FILE BELOW:
# raw_transcripts.gff - the spliced alignments converted into GFF2 format (one transcript per reads)

#clustered algnment coordinates of the reads either as is or after collapsing based on the 5' and 3' coordinates assuming degradation effect 
# clustered_transcripts.gff - The transcripts resulting from the clustering process by cluster_gff.
# clustered_transcripts_collapsed.gff - The transcripts resulting from the clustering with the likely degradation artifacts filtered out.

#clustered algnment coordinates of the polisehed reads either as is or after collapsing based on the 5' and 3' coordinates assuming degradation effect. The grouping of the reads with each other was used to create a consensus sequence devoid of potential basecalling errors. The higher the per transcript coverage the higher the accuracy. The fast sequence of the polished transcripts is provided.   
# polished_transcripts.fas - The fasta sequences of the polished transcripts (one per cluster) produced by polish_clusters.
# polished_transcripts.gff - The alignments of the polished transcripts converted into GFF2 format.
# WE WILL USE ON IGV THE FILE BELOW:
# polished_transcripts_collapsed.gff - The polished transcripts GFF with the likely degradation artifacts filtered out.

#After correcting the reads with each other, the reads were alighned on the genome so that lower coverage reads are furter corrected with the genomic sequence.
# corrected_transcriptome_polished_collapsed.fas - The reference corrected transcriptome generated from the input genome and polished_transcripts_collapsed.gff.

cd ..

# go back to the folder "sequenced_files" to create some alignment statistics of the bam file
cd sequenced_files
#nanoplot of the bam_file
NanoPlot  -o nanoplot_of_summary_bam_file -t 2 --bam /home/wtac/Desktop/new_approach/2019_scripts_final/pipeline_pinfish_analysis_folder/Workspaces/pipeline-pinfish-analysis/alignments/reads_aln_sorted.bam  --maxlength 40000 --plots hex dot
cd ..

#you can also run another alignment statistics software AlignQC
#pip2 install AlignQC
mkdir Align_QC_folder
cd Align_QC_folder
#this you need to run with python 2
#change the alias python=python3 to python=python2 in nano  "~/.bashrc"

#if enough memory you can introduce as the parameter the genome to get more graphs
#alignqc analyze /home/wtac/Desktop/new_approach/2019_scripts_final/pipeline_pinfish_analysis_folder/Workspaces/pipeline-pinfish-analysis/alignments/reads_aln_sorted.bam --threads 1 --specific_tempdir ./tmp -g /home/wtac/Desktop/new_approach/2019_scripts_final/chr_19_mm10_genome_ens75_with_SIRV_and_ERCC/mm10_chr19_ens75_SIRV_original_ERCC.fa --no_transcriptome -o long_reads.alignqc.xhtml --output_folder alignQC.ouput > alignqc.stdout 2>&1

#given memory limitation we can us ethe following
alignqc analyze /home/wtac/Desktop/new_approach/2019_scripts_final/pipeline_pinfish_analysis_folder/Workspaces/pipeline-pinfish-analysis/alignments/reads_aln_sorted.bam --no_genome --threads 1 --specific_tempdir ./tmp --no_transcriptome -o long_reads.alignqc.xhtml --output_folder alignQC.ouput > alignqc.stdout 2>&1
cd ..

#you can also install scanti to analyze the denovo reads
#git clone https://bitbucket.org/ConesaLab/sqanti.git
#cd sqanti/
#sudo chmod a+x *
#pip2 install psutil
#sudo ln -s libudev.so /lib/x86_64-linux-gnu/libudev.so.0
#sudo ln -s libpng16.so /usr/lib/x86_64-linux-gnu/libpng12.so.0

PATH="echo $PATH:/home/wtac/Desktop/new_approach/2019_scripts_final/sqanti/"
mkdir sqanti_analysis_output
cd sqanti_analysis_output
sqanti_qc.py -t 1 -o sqanti_output_analysis --gtf /home/wtac/Desktop/new_approach/2019_scripts_final/pipeline_pinfish_analysis_folder/Workspaces/pipeline-pinfish-analysis/results/polished_transcripts_collapsed.gff /home/wtac/Desktop/new_approach/2019_scripts_final/chr_19_mm10_genome_ens75_with_SIRV_and_ERCC/mm10_chr19_ens75_SIRV_original_ERCC.gtf /home/wtac/Desktop/new_approach/2019_scripts_final/chr_19_mm10_genome_ens75_with_SIRV_and_ERCC/mm10_chr19_ens75_SIRV_original_ERCC.fa
cd ..

#Asses the output on Igv
#write the command  "./../IGV_Linux_2.5.3/igv" on the terminal and press enter
#then:
#1)In the IGV viewer, on the "Genomes" tab select "Load Genome from File" .Then go through the navigator into the folder "chr_19_mm10_genome_ens75_with_SIRV_and_ERCC" that has the fasta sequence of the genome where we aligned to.Select the file "mm10_chr19_ens75_SIRV_original_ERCC.fa" .
#2)Load the bam file mentioned on step (2). On the "File" tab select "Load from File". Then go through the navigator into the folder "pipeline-pinfish-analysis" created into the "alignments" folder.Select the file "reads_aln_sorted.bam" .
#3)Load the gtf file with the gene models. On the "File" tab select "Load from File". Then go through the navigator into the folder "chr_19_mm10_genome_ens75_with_SIRV_and_ERCC".Select the file "mm10_chr19_ens75_SIRV_original_ERCC.gtf".
#5)Load the de novo isoform models from the collapsed nanopore reads from the  file "polished_transcripts_collapsed.gff" inside the folder "results" in the folder "pipeline-pinfish-analysis". 
#6) In the "View" tab, under the "Preferences" option go to the "Alignments" tabs and check whether the box next to "Downsample reads" is unchecked
#7) On the IGV viewer on the box next to "Go" word put genomic coordinates for example "19:29,358,128-29,363,269" to go to the genomic region of interest.


#compare the new models with the known ones
#download from: http://ccb.jhu.edu/software/stringtie/dl/
#

cd 
PATH="echo $PATH:/home/wtac/Desktop/new_approach/2019_scripts_final/gffcompare-0.11.2.Linux_x86_64/"
gffcompare -r /home/wtac/Desktop/new_approach/2019_scripts_final/chr_19_mm10_genome_ens75_with_SIRV_and_ERCC/mm10_chr19_ens75_SIRV_original_ERCC.gtf /home/wtac/Desktop/new_approach/2019_scripts_final/pipeline_pinfish_analysis_folder/Workspaces/pipeline-pinfish-analysis/results/polished_transcripts_collapsed.gff


gffcompare -R -Q -o FOR_PRECISION_STATS -r /home/wtac/Desktop/new_approach/2019_scripts_final/chr_19_mm10_genome_ens75_with_SIRV_and_ERCC/mm10_chr19_ens75_SIRV_original_ERCC.gtf /home/wtac/Desktop/new_approach/2019_scripts_final/pipeline_pinfish_analysis_folder/Workspaces/pipeline-pinfish-analysis/results/polished_transcripts_collapsed.gff


cd ..

#================================================================================================================

#In order to calculate the gene and isoform abundance we will use the "pipeline-transcriptome-de" from nanoporetech    

## General pipeline parameters:

# pipeline: "pipeline-transcriptome-de_phe" -->Name of the pipeline
# workdir_top: "Workspaces" -->ABSOLUTE path to directory holding the working directory
# resdir: "results" -->Results directory
# repo: "https://github.com/nanoporetech/pipeline-transcriptome-de" -->Repository URL

## Pipeline-specific parameters:
# transcriptome: "/home/wtac/Desktop/new_approach/s1_WT/chromosomes/Mus_musculus.GRCm38.cdna.all.fa" -->Transcriptome fasta download from the Ensembl website (ftp://ftp.ensembl.org/pub/release-96/fasta/mus_musculus/cdna/)
# annotation: "/home/wtac/Desktop/new_approach/s1_WT/chromosomes/Mus_musculus.GRCm38.96.gtf" -->--> Annotation GFF/GTF file downloaded from the Ensembl website (ftp://ftp.ensembl.org/pub/release-96/gtf/mus_musculus/)

# The list of control samples as many as available:
#C1: "/home/wtac/Desktop/new_approach/multiple_files/10000_s1_WT_pass_reads.fastq"
#C2: "/home/wtac/Desktop/new_approach/multiple_files/10000_s2_WT_pass_reads.fastq"

# The list of treated samples as many as available:
#IR1: "/home/wtac/Desktop/new_approach/multiple_files/10000_s4_scramble_pass_reads.fastq"
#IR2: "/home/wtac/Desktop/new_approach/multiple_files/10000_s5_scramble_pass_reads.fastq"

# minimap_index_opts: "" -->Minimap2 indexing options. Leave it as default (empty)
# minimap2_opts: "" -->Minimap2 mapping options . Leave it as default (empty)

# maximum_secondary: 100 -->Maximum secondary alignments. Because we align on the transcriptome there can be multiple similar isoforms where the 
# secondary_score_ratio: 1.0 -->Secondary score ratio (-p for minimap2). Retain up to the above number top secondary mappings if their chaining scores are higher than -p [=1] of their corresponding primary mappings. Practically given that we used p=1 tghis means that we accept secondary alignments thata re as good as the primary alignment.


# salmon_libtype: "U" -->Salmon library type. This is an option of the Salmon aligner and specifies whether the protocol is stranded (parameter value:"S") or unstranded (parameter value:"U"); As we have transfrom all the mRNA in their cDNA equivalent we will select the parameter: "S". 

# Count filtering options - customize these according to your experimental design:

# min_samps_gene_expr: 1 -->Genes expressed in minimum this many samples
# min_samps_feature_expr: 1 -->Transcripts expressed in minimum this many samples
# min_gene_expr: 1 -->Minimum gene counts
# min_feature_expr: 1 -->Minimum transcript counts

# threads: 50 --> number of threads 


#The software outputs the following files in the "Workspaces" folder:
#    alignments/*.bam - unsorted transcriptome alignments (input to salmon).
#    alignments_sorted/*.bam - sorted and indexed transcriptome alignments.
#    counts - counts generated by salmon.
#    merged/all_counts.tsv - the transcript count table including all samples.
#    merged/all_counts_filtered.tsv - the transcript count table including all samples after filtering.
#    merged//all_gene_counts.tsv - the gene count table including all samples.
#    de_analysis/coldata.tsv - the condition table used to build model matrix.
#    de_analysis/de_params.tsv - analysis parameters generated from config.yml.
#    de_analysis/results_dge.tsv and de_analysis/results_dge.pdf- results of edgeR differential gene expression analysis.
#    de_analysis/results_dtu_gene.tsv, de_analysis/results_dtu_transcript.tsv and de_analysis/results_dtu.pdf - results of differential transcript usage by DEXSeq.
#    de_analysis/results_dtu_stageR.tsv - results of the stageR analysis of the DEXSeq output.
#    de_analysis/dtu_plots.pdf - DTU results plot based on the stageR results and filtered counts.


#IMPORTANT:We note here that the software will produce the count files of the transcriptomes and teh genes only if  the differential expression analysis is succesful which indicates that it needs replicates in both the control and the sample files otherwise the variance will not be calculated properly and the software will crush. If only 1 file will be used maybe dummy files can be intrduced that will permit the sorftware to run into completion.   

#to address the lack of other samples we will subsample reads from the bamfile. This will  create some artificial variance that will be necessary to complete the execution of the software. 
#-s is the random seed

cd sequenced_files
seqtk sample -s100 adapter_removed_full_length_above_q8_pass_reads.fastq 200000 >dummy_sample1.fastq
seqtk sample -s95 adapter_removed_full_length_above_q8_pass_reads.fastq 150000 >dummy_sample2.fastq
seqtk sample -s90 adapter_removed_full_length_above_q8_pass_reads.fastq 100000 >dummy_sample3.fastq
cd ..

#Then you can run the file of interest

#head -5000 quality_8_filtered_pass_reads.fastq >5000_quality_8_filtered_pass_reads.fast
cd pipeline-transcriptome-de
snakemake --use-conda -j 1 all
cd ..

#=======================================

