WTAC_Transcriptomics_2019

cDNA library preparationa and sequencing with Oxford Nanopore Technologies














Basecalling with Guppy

Run Guppy on Linux

Launch a terminal on your system. For example, when using Ubuntu, you can type Ctrl + Alt + T.

You will run all Guppy operations from here:

guppy_basecaller --input_path … --save_path … --flowcell … --kit -q 0 –num_callers 4 –qscore-filtering –min-qscore 7 -r
