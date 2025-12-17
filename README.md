
# numap R package

**Maintainer**: Ji-Ping Wang, <jzwang@northwestern.edu>

**Reference**: 

- Li, K., Tu, Y., Voong, L.N., Lu, Y., Kendall, B., Ma, X., Pui, S.L., Tao, M., Wang, J.P., and Wang, X., Differential nucleosome organization in human interphase and metaphase chromosomes, 2025
- Brogaard, K., Xi, L., Wang, J.-P., and Widom, J., A Map of Nucleosome Positions in Yeast at base-pair resolution, Nature, 2012

## What is `numap`?

`numap` is an R pipeline designed to process chemical map of nucleosome positioning data. It provides a range of functions from processing BAM files, computing the unique and redundant nucleosome maps, and nucleosome occupancy. The following are the major functions of this package.

---

## NCP_wrapper1

A wrapper function to compute the pipeline including parsing BAM, single-template deconvolution, unique and redundant map calling, and occupancy calculation.

### Usage

```r
NCP_wrapper1(
  bam_file,           # name of a single BAM file
  paired=TRUE,        # paired or single-end
  result_path,        # folder where results are saved
  chromosome_path,    # folder where chromosome files are
  temp1=NULL,         # single-template file
  center=TRUE,        # center or uniformly weighted occupancy
  unique_map_wsize=120, # window size for unique map
  wsize=73,           # nucleosome occupancy window size
  parse_skip=FALSE,   # skip parsing BAM file if TRUE
  chrom_to_include=NULL # vector of chromosome names to include
)
```

### Example

```r
bam_file1 <- "/Users/jon/Rpack/cheMNuP_test/genome/bam/yeast_77_78/yeast_77_78.bam"
result_path1 <- "/Users/jon/Rpack/cheMNuP_test/genome/bam/yeast_77_78"
chromosome_path1 <- "/Users/jon/Rpack/cheMNuP_test/genome/bam/"

NCP_wrapper1(
  bam_file = bam_file1,
  paired = TRUE,
  temp1 = NULL,
  chromosome_path = chromosome_path1,
  result_path = result_path1,
  center = FALSE,
  parse_skip = FALSE
)
```

---

## NCP_wrapper4

A wrapper function for 4-template deconvolution, unique and redundant map calling, and occupancy calculation.

### Usage

```r
NCP_wrapper4(
  bam_file,           # name of a single BAM file
  paired=TRUE,        # paired or single-end
  result_path,        # folder where results are saved
  chromosome_path,    # folder where chromosome files are
  temp4=NULL,         # four-template file
  center=TRUE,        # center or uniformly weighted occupancy
  unique_map_wsize=120, # window size for unique map
  wsize=73,           # nucleosome occupancy window size
  parse_skip=FALSE,   # skip parsing BAM file if TRUE
  chrom_to_include=NULL # vector of chromosome names to include
)
```

### Example

```r
NCP_wrapper4(
  bam_file = bam_file1,
  paired = TRUE,
  temp4 = NULL,
  chromosome_path = chromosome_path1,
  result_path = result_path1,
  center = FALSE,
  parse_skip = FALSE
)
```

---

## NCP1

Main function to perform deconvolution algorithm using 1 template.

### Usage

```r
NCP1(
  result_path,         # directory where results are stored
  chromosome_path,     # path where chromosomes are stored, one file per chromosome
  temp1=NULL,          # file containing the template
  unique_map_wsize=120 # distance to define local maximum for unique map
)
```

---

## NCP4

Main function to perform deconvolution algorithm using 4 templates.

### Usage

```r
NCP4(
  result_path,         # directory where results are stored
  chromosome_path,     # path where chromosomes are stored, one file per chromosome
  temp4=NULL,          # file containing the 4-template
  unique_map_wsize=120 # distance to define local maximum for unique map
)
```

---

## occupancy_Mnase_paired_par

Compute nucleosome occupancy based on MNase maps in parallel.

### Usage

```r
occupancy_Mnase_paired_par(
  reads_file_list,       # list of parsed BAM files from parse_bam_paired_Mnase
  result_path,           # directory to save results
  chrom_length=NULL,     # vector of chromosome lengths, if NULL use max end + wsize
  center=TRUE,           # calculate center-weighted or uniformly weighted occupancy
  wsize=73,              # window size for weighted occupancy, ignored if center=FALSE
  insert_length_low,     # lower limit of read length to include
  insert_length_high,    # upper limit of read length to include
  cores=NULL,            # number of cores to use in parallel
  maxhits=Inf,           # maximum hits a read can have to be included
  nameExt=NULL           # extension to add to output file
)
```

---

## occupancy_Mnase_paired

Compute nucleosome occupancy based on MNase maps.

### Usage

```r
occupancy_Mnase_paired(
  reads_file_list,       # list of parsed BAM files from parse_bam_paired_Mnase
  result_path,           # directory to save results
  chrom_length=NULL,     # vector of chromosome lengths, if NULL use max end + wsize
  center=TRUE,           # calculate center-weighted or uniformly weighted occupancy
  wsize=73,              # window size for weighted occupancy, ignored if center=FALSE
  insert_length_low,     # lower limit of read length to include
  insert_length_high,    # upper limit of read length to include
  maxhits=Inf,           # maximum hits a read can have to be included
  nameExt=NULL           # extension to add to output file
)
```

---

## occupancy

Compute nucleosome occupancy.

### Usage

```r
occupancy(
  redundant_map_list,  # list of redundant maps by chromosome
  result_path,         # directory to save results
  chrom_length=NULL,   # length of each chromosome, if NULL use max coordinate + wsize
  T4=FALSE,            # FALSE for single-template, TRUE for 4-template
  center=TRUE,         # center-weighted or uniformly weighted occupancy
  wsize=73             # range of weights, +/- wsize
)
```

---

## parse_bam_paired_Mnase

Parse BAM file of paired MNase map.

### Usage

```r
parse_bam_paired_Mnase(
  bam_file,          # BAM file from paired-end alignment
  result_path,       # directory to save parsed BAM
  chrom_to_include   # vector of chromosome names to process
)
```

---

## parse_bam_paired

Parse paired-end chemical mapping BAM file.

### Usage

```r
parse_bam_paired(
  bam_file,          # BAM file from paired-end reads
  result_path,       # directory to save results
  chrom_to_include   # vector of chromosome names to process
)
```

---

## parse_bam_single

Parse single-end chemical mapping BAM file.

### Usage

```r
parse_bam_single(
  bam_file,          # single-end BAM file
  result_path,       # directory to save results
  chrom_to_include   # vector of chromosome names to process
)
```

---

## plotAATT

Plot AA/TT/TA/AT frequency based on a nucleosome map.

### Usage

```r
plotAATT(
  chromosome_path,   # path to chromosome FASTA files
  center_list,       # list of nucleosome maps by chromosome
  result_path,       # directory to save results
  NCP_thresh=0,      # threshold to select nucleosomes from unique map
  wsize=73           # window size to plot AA/TT/TA/AT frequency
)
```

---

## TSS_occu

Calculate occupancy around TSS.

### Usage

```r
TSS_occu(
  TSS_file,    # tab- or space-delimited file with chr, strand, start (1=Watson, 2=Crick)
  occu_list,   # occupancy file output from occupancy function
  wsize=1000   # window size, default +/- 1000 bp
)
```
