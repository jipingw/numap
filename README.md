numap R package
================

**Maintainer**: Ji-Ping Wang, \<<jzwang@northwestern.edu>\>

**Reference**: 

-Li, K., Tu, Y., Voong, L.N., Lu, Y., Kendall, B., Ma, X., Pui, S.L., Tao, M., Wang, J.P., and  Wang, X.,  Differential nucleosome organization in human interphase and metaphase chromosomes, 2025

-Brogaard, K., Xi, L., Wang, J.-P., and Widom, J., A Map of Nucleosome Positions in Yeast at base-pair resolution, Nature, 2012

## What is `numap`?
`numap` is an R pipeline designed to process chemical map of nucleosome positioning data. It provides a range of functions from processing of bam file, compute the unique and redundant nucleosome maps and nucleosome occupancy. The following are the major functions of this package.


### NCP_wrapper1

A wrapper function to compute the pipeline including parseing bam, single-template deconvolution, unique and redundant map calling and occupancy calculation.

#### Usage

``` R
NCP_wrapper1(bam_file,paired=TRUE, result_path, chromosome_path,
                      temp1=NULL,center=TRUE,unique_map_wsize=120, 
wsize=73,parse_skip=FALSE,chrom_to_include =
                 NULL)
```

#### Arguments

  --------------------------------------- ----------------------------------------------------------------------
  `bam_file`{#bam_file}                   name of a single bam file
  `paired`{#paired}                       paired or single-end
  `result_path`{#result_path}             the folder where the results files are saved
  `chromosome_path`{#chromosome_path}     the folder where the chromosome files are
  `temp1`{#temp1}                         single template file
  `center`{#center}                       center or uniformly weighted occupancy to be calculate
  `unique_map_wsize`{#unique_map_wsize}   the window size to define the local maximum for the unique map
  `wsize`{#wsize}                         nucleosome occupancy window size for weights
  `parse_skip`{#parse_skip}               skip parse bam file if TRUE
  `chrom_to_include`{#chrom_to_include}   a vector listing the chromosome names of which the data to be output
  --------------------------------------- ----------------------------------------------------------------------

#### Examples

``` {r eval=FALSE}
bam_file1="/Users/jon/Rpack/cheMNuP_test/genome/bam/yeast_77_78/yeast_77_78.bam"
result_path1="/Users/jon/Rpack/cheMNuP_test/genome/bam/yeast_77_78"
temp1="/Users/jon/Dropbox/chemHuman/scripts/other/template1.txt"
temp4="/Users/jon/Dropbox/chemHuman/scripts/other/template4.txt"
chromosome_path1="/Users/jon/Rpack/cheMNuP_test/genome/bam/"
#sourceCpp("/Users/jon/Rpack/numap/src/estNCP_Rcpp.cpp")

NCP_wrapper1(bam_file=bam_file1,paired=TRUE,temp1=NULL,chromosome_path=chromosome_path1,result_path=result_path1,center=FALSE,parse_skip=FALSE)
```

---

### NCP_wrapper4
A wrapper function to compute the pipeline including parseing bam, 4-template deconvolution, unique and redundant map calling and occupancy calculation.


#### Usage

``` {r eval=FALSE}
NCP_wrapper4(bam_file,paired=TRUE, result_path, chromosome_path,
                      temp4=NULL,center=TRUE,unique_map_wsize=120, 
wsize=73,parse_skip=FALSE,chrom_to_include =
                 NULL)
```

#### Arguments

  --------------------------------------- -----------------------------------------------------------------------
  `bam_file`{#bam_file}                   name of a single bam file
  `paired`{#paired}                       paired or single-end
  `result_path`{#result_path}             the folder where the results files are saved
  `chromosome_path`{#chromosome_path}     the folder where the chromosome files are
  `temp4`{#temp4}                         four template file
  `center`{#center}                       center or uniformly weighted occupancy to be calculate
  `unique_map_wsize`{#unique_map_wsize}   the window size to define the local maximum for the unique map
  `wsize`{#wsize}                         nucleosome occupancy window size for weights
  `parse_skip`{#parse_skip}               skip parse bam file if TRUE
  `chrom_to_include`{#chrom_to_include}   a vector listing chromosome names of which the mapping data to output
  --------------------------------------- -----------------------------------------------------------------------

#### Examples

```{r eval=FALSE}
bam_file1="/Users/jon/Rpack/cheMNuP_test/genome/bam/yeast_77_78/yeast_77_78.bam"
result_path1="/Users/jon/Rpack/cheMNuP_test/genome/bam/yeast_77_78"
temp1="/Users/jon/Dropbox/chemHuman/scripts/other/template1.txt"
temp4="/Users/jon/Dropbox/chemHuman/scripts/other/template4.txt"
chromosome_path1="/Users/jon/Rpack/cheMNuP_test/genome/bam/"
#sourceCpp("/Users/jon/Rpack/numap/src/estNCP_Rcpp.cpp")

NCP_wrapper4(bam_file=bam_file1,paired=TRUE,temp4=NULL,chromosome_path=chromosome_path1,result_path=result_path1,center=FALSE,parse_skip=FALSE)

```

---

### NCP1
Main function to perform deconvolution algorithm using 4 templates.

#### Usage

``` R
NCP1(result_path, chromosome_path, temp1=NULL,unique_map_wsize=120)
```

#### Arguments

  --------------------------------------- ------------------------------------------------------------------------
  `result_path`{#result_path}             directory where the results are stored
  `chromosome_path`{#chromosome_path}     path to where the chromosomes are stored, one file for each chromosome
  `temp1`{#temp1}                         file contains the template
  `unique_map_wsize`{#unique_map_wsize}   distance to define the local maximum for unique map
  --------------------------------------- ------------------------------------------------------------------------

#### Value

NCP1 generates the NCP scores and the unique and redudnant maps


---

### NCP4

Main function to perform deconvolution algorithm using 4 templates.

#### Usage

``` R
NCP4(result_path, chromosome_path, temp4=NULL,unique_map_wsize=120)
```

#### Arguments

  --------------------------------------- ------------------------------------------------------------------------
  `result_path`{#result_path}             directory where the results are stored
  `chromosome_path`{#chromosome_path}     path to where the chromosomes are stored, one file for each chromosome
  `temp4`{#temp4}                         file contains the 4-template
  `unique_map_wsize`{#unique_map_wsize}   distance to define the local maximum for unique map
  --------------------------------------- ------------------------------------------------------------------------

#### Value

NCP4 generates the NCP scores and the unique and redundant maps

#### occupancy_Mnase_paired_par

This function compute the nucleosome occupancy based on MNase maps in
parallel.

#### Usage

``` R
occupancy_Mnase_paired_par(reads_file_list, result_path, chrom_length=NULL,center=TRUE,wsize=73, insert_length_low, insert_length_high,cores=NULL, maxhits=Inf,nameExt=NULL)
```

#### Arguments

  ------------------------------------------- -----------------------------------------------------------------------------------------------------
  `reads_file_list`{#reads_file_list}         list of parsed bam file from \'parse_bam_paired_Mnase function
  `result_path`{#result_path}                 directory to save the results
  `chrom_length`{#chrom_length}               vector to contain the chromosome length. if \'NULL\', the maximum end position + wsize will be used
  `center`{#center}                           whether to calcualte center weighted or uniformly weighted occupancy
  `wsize`{#wsize}                             the widown size to calcualte the center weighted occupany. Ignored if \'center=FALSE\'
  `insert_length_low`{#insert_length_low}     the lower limit of read length to keep for calculation of occupancy
  `insert_length_high`{#insert_length_high}   the upper limit of read length to keep for calculation of occupancy
  `cores`{#cores}                             the number of cores to use in parallel
  `maxhits`{#maxhits}                         the maximum of hits that a read can have to be included in occupancy calculation
  `nameExt`{#nameExt}                         the extension to be added to the occupancy output file name
  ------------------------------------------- -----------------------------------------------------------------------------------------------------

#### Value

Output occupancy value by chromosomes.

---

### occupancy_Mnase_paired

This function compute the nucleosome occupancy based on MNase maps.

#### Usage

``` R
occupancy_Mnase_paired(reads_file_list, result_path, chrom_length=NULL, center=TRUE,wsize=73, insert_length_low, insert_length_high,maxhits=Inf,nameExt=NULL)
```

#### Arguments

  ------------------------------------------- -----------------------------------------------------------------------------------------------------
  `reads_file_list`{#reads_file_list}         list of parsed bam file from \'parse_bam_paired_Mnase function
  `result_path`{#result_path}                 directory to save the results
  `chrom_length`{#chrom_length}               vector to contain the chromosome length. if \'NULL\', the maximum end position + wsize will be used
  `center`{#center}                           whether to calcualte center weighted or uniformly weighted occupancy
  `wsize`{#wsize}                             the widown size to calcualte the center weighted occupany. Ignored if \'center=FALSE\'
  `insert_length_low`{#insert_length_low}     the lower limit of read length to keep for calculation of occupancy
  `insert_length_high`{#insert_length_high}   the upper limit of read length to keep for calculation of occupancy
  `maxhits`{#maxhits}                         the maximum of hits that a read can have to be included in occupancy calculation
  `nameExt`{#nameExt}                         the extension to be added to the occupancy output file name
  ------------------------------------------- -----------------------------------------------------------------------------------------------------

#### Value

Output occupancy value by chromosomes.

---

### occupancy

function to compute nucleosome occupancy

#### Usage

``` R
occupancy(redundant_map_list,result_path, chrom_length=NULL,T4=FALSE,
center=TRUE,wsize=73)
```

#### Arguments

  ------------------------------------------- --------------------------------------------------------------------------------------------------------------------------------------------
  `redundant_map_list`{#redundant_map_list}   a list of redundant maps by chromosomes
  `result_path`{#result_path}                 directory where the results to be saved
  `chrom_length`{#chrom_length}               the length of each chromosomes. If not provided, the occupancy is calculated based on the maximum coordinate of the redundant map + wsize.
  `T4`{#T4}                                   T4=false to compute single-template, true for 4-template
  `center`{#center}                           center-weighted or uniformly weighted occupancy
  `wsize`{#wsize}                             range of weights, +/- wsize
  ------------------------------------------- --------------------------------------------------------------------------------------------------------------------------------------------

#### Value

This function outputs calculated occupancy per chromosome in files.

---

### parse_bam_paired_Mnase

This function parses bam file of paired MNase map.

#### Usage

``` R
parse_bam_paired_Mnase(bam_file, result_path,chrom_to_include)
```

#### Arguments

  --------------------------------------- -----------------------------------------------------------------------------
  `bam_file`{#bam_file}                   the bam file from paired end alignment
  `result_path`{#result_path}             directory path to save the parsed bam file
  `chrom_to_include`{#chrom_to_include}   a vector listing chromosome names of which the mapping data to be processed
  --------------------------------------- -----------------------------------------------------------------------------

#### Value

output parsed bam file, each row contains start and end positions.

---

### parse_bam_paired
function to parse paired-end chemical mapping bam file

#### Usage

``` R
parse_bam_paired(bam_file,result_path,chrom_to_include)
```

#### Arguments

  --------------------------------------- -------------------------------------------------------------------
  `bam_file`{#bam_file}                   bam file from paired-end reads
  `result_path`{#result_path}             directory where results are stored
  `chrom_to_include`{#chrom_to_include}   a vector listing the chromosome names of which the data to output
  --------------------------------------- -------------------------------------------------------------------

#### Value

output cleavage frequency at each genomic location on watson and crick
strands separately for each chromosome.

---

### parse_bam_single
function to parse single-end chemical mapping bam file

#### Usage

``` R
parse_bam_single(bam_file,result_path,chrom_to_include)
```

#### Arguments

  --------------------------------------- -------------------------------------------------------------------
  `bam_file`{#bam_file}                   single-end paired end bam file
  `result_path`{#result_path}             directory where results are stored
  `chrom_to_include`{#chrom_to_include}   a vector listing the chromosome names of which the data to output
  --------------------------------------- -------------------------------------------------------------------

#### Value

output cleavage frequency at each genomic location on watson and crick
strands separately for each chromosome.
---

### plotAATT
Function to plot AA/TT/TA/AT fequency based on a nucleosome map

#### Usage

``` R
plotAATT(chromosome_path,center_list,result_path,NCP_thresh=0, wsize=73)
```

#### Arguments

  ------------------------------------- ----------------------------------------------------------------------
  `chromosome_path`{#chromosome_path}   path to chromosome fasta files
  `result_path`{#result_path}           directory where the results to be saved
  `center_list`{#center_list}           list of nucleosome maps by chromosomes
  `NCP_thresh`{#NCP_thresh}             NCP_thresh to define which nucleosome to be used from the unique map
  `wsize`{#wsize}                       window size to plot the AATT frequency
  ------------------------------------- ----------------------------------------------------------------------

#### Value

It returns AA/TT/TA/AT frequency plot in-screen around nucleosome
centers. The function also returns a vector of AA/TT/TA/AT frequency in
the defined region.
---
### TSS_occu
This function calculates occupancy around TSS.

#### Usage

``` R
TSS_occu(TSS_file, occu_list, wsize=1000)
```

#### Arguments

  ------------------------- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
  `TSS_file`{#TSS_file}     TSS_file must tab or space_delimited file and contain three columns, chr, strand and start (start is TSS start), 1 for watson and 2 for crick strand. chr must match the chromosome names in Occupancy file.
  `occu_list`{#occu_list}   occupancy file outputted by \'occupancy\' function
  `wsize`{#wsize}           window size. default is +/-100 bp.
  ------------------------- --------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

