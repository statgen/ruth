# RUTH - Robust Unified Hardy-Weinberg Equilibrium Test

`ruth` is a software to perform robust unified Hardy-Weinberg Equilbrium (HWE) tests for sequence-based genotypes under population structure.

### Quick Overview

Inputs: 
 * A genotype file in VCF or BCF format, with either genotype likelihoods (PLs or GLs) or best-guess genotypes (GTs)
 * Principal components (PCs) or some other ancestry summary statistics for the samples in the VCF or BCF file
 
Outputs:
 * Robust Hardy-Weinberg Equilibrium statistics, which accounts for the effects of population stratification by using principal components

### Introduction

#### Overview

'ruth' (Robust Unified Test for HWE) uses information from genotypes and principal components to perform either a likelihood ratio test or a score test to estimate variants' deviation from HWE after adjusting for population structure. 

### Tips for running

 * The user needs to have a **genotype file in VCF or BCF format** and **estimated PCs** for the samples 
     * If PCs are not available, they can be calculated if you have access to the aligned sequences (BAM or CRAM files) using VerifyBamID2, available here: https://github.com/Griffan/VerifyBamID
 * To decrease the size of the output file, use the **--site-only** option to suppress outputting individual-level genotypes
 * If available, we recommend using genotype likelihoods (either "--field PL" or "--field GL")
 * We currently recommend setting lambda to 0 (**--lambda 0**), and using the likelihood ratio EM test (**--lrt-em**)

### Installing ruth

<pre>
$ mkdir build

$ cd build

$ cmake ..
</pre>

In case any required libraries is missing, you may specify customized installing path by replacing "cmake .." with:

<pre>
For libhts:
  - $ cmake -DHTS_INCLUDE_DIRS=/hts_absolute_path/ -DHTS_LIBRARIES=/hts_absolute_path/libhts.a ..

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
All software use a self-documentation utility. You can run each utility with -man or -help option to see the command line usages. Also, we offer some general practice with a tutorial example (data available here: [TBA])

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

</pre>

### Interpretation of output files

 * The definitions of the added INFO fields can be found in the header of the output VCF or BCF file
 * The statistic of interest is **HWE_SLP_I**, which is the signed log10 P-value of the HWE test with individual-specific allele frequencies, adjusted for population structure
    * HWE_SLP_I < 0 indicates **an excess of heterozygotes**
    * HWE_SLP_I > 0 indicates **heterozygote depletion**
    * An excess of heterozygotes can be a telltale sign of certain types of technical artefacts
 * Any P-value threshold represents a tradeoff between sensitivity and specificity
    * Using a more stringent threshold will decrease false positives but increase false negatives
    * Using a less stringent threshold will have the opposite effect
 * With low coverage data, a slightly more stringent threshold (for example, P < 1e-4) can help with reducing false positives
 * With high coverage data, a less stringent threshold (for example, P < 0.01 or P < 0.001) can lead to improved power while maintaining good false positive performance
  
