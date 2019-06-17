# RNA Transcriptomics
18-29th June 2019, Wellcome Genome Campus, Hinxton, Cambridge, CB10 1RQ

![alt text](https://github.com/abayega/Courses-and-Practicals/blob/master/WTAC_2019/images/mcgill%20logo.png)

## cDNA library preparation and sequencing with Oxford Nanopore Technologies



## Basecalling with Guppy

![Short Tutorial](https://github.com/abayega/Courses-and-Practicals/blob/master/WTAC_2019/Summary%20on%20ONT%20Basecalling.pdf)

#### Run Guppy on Linux

#### Launch a terminal on your system. For example, when using Ubuntu, you can type Ctrl + Alt + T.

### You will run all Guppy operations from here:
```
guppy_basecaller --input_path … --save_path … --flowcell … --kit -q 0 –num_callers 4 –qscore-filtering –min-qscore 7 -r
```

## Reads quality control checks
Initially we will start with the sequencing run statistics derived from the summary file or from the fastq file
go inside the "summary file" folder and process the "sequencing_summary.txt" file
```
cd summary_file_folder
```
Run NanoPlot to derive sequencing summary stats
```
NanoPlot -o nanoplot_of_summary_file --readtype 1D --summary sequencing_summary.txt --loglength
```
```
cd ..
```

Here we will use Primer-chop to identify full length reads based on the presence or absence of the adaptor sequences in the beggining or the end of the reads.
```
cd sequenced_files
```
```
./../primer-chop/primer-chop -P2 ../primers.fa above_q8_pass_reads.fastq primer_chop_output
```
Here we will try to find how many reads have been reported in each one of the produced files.
```
wc -l primer_chop_output/good-fwd.fa
wc -l primer_chop_output/good-rev.fa
wc -l primer_chop_output/no_head-fwd.fa
wc -l primer_chop_output/no_head-rev.fa
wc -l primer_chop_output/no_tail-fwd.fa
wc -l primer_chop_output/no_tail-rev.fa
```

take the files that contain  either the  sence or antisense reads that have adaptors present at both end of the molecules
```
cat primer_chop_output/good-fwd.fa primer_chop_output/good-rev.fa >adapter_removed_full_length_above_q8_pass_reads.fa
```

## De novo genome guided transcriptome assembly
```
PATH="$PATH:/home/wtac/Desktop/new_approach/BBMap_38.50b/bbmap/"
reformat.sh in=adapter_removed_full_length_above_q8_pass_reads.fa out=adapter_removed_full_length_above_q8_pass_reads.fastq qfake=12
cd ..
```
we align on the genome in order to identify new spliced isoforms
add minimap2 in path!!!
```
PATH="$PATH:/home/wtac/Desktop/new_approach/s1_WT/minimap2-2.17_x64-linux"
cd pipeline_pinfish_analysis_folder #GO INTO THIS FOLDER
snakemake --use-conda -j 1 all
```

## Gene expression counts

