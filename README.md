# ruth - Robust Unified Hardy-Weinberg Equilibrium Test

`ruth` is a software to perform robust unified Hardy-Weinberg Equilbrium (HWE) test for sequence-based genotypes under population structure.

### Quick Overview

TBA

### Introduction

#### Overview

TBA

### Tips for running

TBA

### Installing ruth

<pre>
$ mkdir build

$ cd build

$ cmake ..
</pre>

In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

<pre>
For libhts:
  - $ cmake -DHTS_INCLUDE_DIRS=/hts_absolute_path/include/  -DHTS_LIBRARIES=/hts_absolute_path/lib/libhts.a ..

For bzip2:
  - $ cmake -DBZIP2_INCLUDE_DIRS=/bzip2_absolute_path/include/ -DBZIP2_LIBRARIES=/bzip2_absolute_path/lib/libbz2.a ..

For lzma:
  - $ cmake -DLZMA_INCLUDE_DIRS=/lzma_absolute_path/include/ -DLZMA_LIBRARIES=/lzma_absolute_path/lib/liblzma.a ..
</pre>

Finally, to build the binary, run

<pre>
$ make
</pre>

### Using ruth
All softwares use a self-documentation utility. You can run each utility with -man or -help option to see the command line usages. Also, we offer some general practice with an example in tutorial (data is available here: [TBA])

<pre>
$(RUTH_HOME)/bin/ruth --vcf [Input VCF file] --evec [Input EigenVector] --out [Output]
</pre>

The detailed usage is also pasted below.

<pre>
Input Options
   --evec        [STR: ]             : (REQUIRED) Name of eigenvector file, where each line contains [SAMPLE_ID] [PC1] [PC2] ..... The number of PCs could be larger than parameters specified by --num-PC
   --vcf         [STR: ]             : (REQUIRED) Input VCF/BCF file
   --thin        [FLT: 1.00]         : Probability to randomly sample variants from BCF
   --seed        [INT: 0]            : Random seed to set (default is to use the clock time)
   --num-pc      [INT: 4]            : Number of principal componentds to be used from the file specified by --evec 
   --field       [STR: ]             : FORMAT field in VCF to extract the genotype likelihood or genotypes from. Only PL, GL, GT are allowed currently
   --gt-error    [FLT: 5.0e-03]      : Error rates for GT field when --field GT option is used. Ignored for other fields
   --lambda      [FLT: 1.00]         : Max lambda parameter

Output Options
   --out         [STR: ]             : (REQUIRED) Output VCF file to write with ISHWEZ and ISIBC statistics and IF format field
   --skip-if     [FLG: OFF]          : Skip writing individual-specific allele frequency for each sample in output VCF/BCF
   --skip-info   [FLG: OFF]          : Skip updating INFO field for each sample in output VCF/BCF
   --site-only   [FLG: OFF]          : Do not write genotype information, and writes only site information (up to INFO field) in output VCF/BCF
   --nelder-mead [FLG: OFF]          : Use Nelder-Mead algorithm (instead of EM) when estimating individual-specific allele frequencies
   --lrt-test    [FLG: OFF]          : Use Likelihood-ratio test with Nelder-Mead algorithm (instead of score test) for performing HWE test
   --lrt-em      [FLG: OFF]          : Use Likelihood-ratio test with EM algorithm (instead of score test) for performing HWE test

Samples to focus on
   --sm-list     [STR: ]             : A file containg the list of sample IDs to subset

Parameters for sex chromosomes
   --sex-map     [STR: ]             : Sex map file in PED format or tsv file with [ID,SEX in X ploidy]
   --x-label     [STR: X]            : Label for X chromosome
   --y-label     [STR: Y]            : Label for Y chromosome
   --mt-label    [STR: MT]           : Label for MT chromosome
   --x-start     [INT: 2699520]      : Start coordinate of non-PAR X region
   --x-stop      [INT: 154931044]    : Stop coordinate of non-PAR X region

Options to specify when chunking is used
   --ref         [STR: ]             : Reference FASTA file name (required only when chunking is used)
   --unit        [INT: 2147483647]   : Chunking unit in bp (specify only with --ref together
   --interval    [STR: ]             : Interval file name used for chunking (specify only when chunking is used without --ref
   --region      [STR: ]             : Target region to focus on

Options for input SAM/BAM/CRAM 
   --sam [STR: ] : Input SAM/BAM/CRAM file. Must be sorted by coordinates and indexed
   --tag-group [STR: CB] : Tag representing readgroup or cell barcodes, in the case to partition the BAM file into multiple groups. For 10x genomics, use CB
   --tag-UMI [STR: UB] : Tag representing UMIs. For 10x genomiucs, use UB

Options for input VCF/BCF
   --vcf [STR: ] : Input VCF/BCF file, containing the AC and AN field
   --sm [V_STR: ] : List of sample IDs to compare to (default: use all)
   --sm-list [STR: ] : File containing the list of sample IDs to compare
</pre>

### Interpretation of output files

TBA
  


