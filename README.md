# Download gene expression data from GEO database
Author: Weiyang Tao

Email: weiyangtao1513@gmail.com

Date: 2023-07-28

------

This package is to help you download gene expression or raw fastq files and meta data (i.e. pheno data) from NCBI Gene Expression Omnibus (GEO) database.

There are two types of gene expression data sets: (1) microarray; (2) RNA-seq.

## Microarray data
For microarray data, GEO database usually stores `Series Matrix File(s)` and sometimes a supplementary file called `GSE*****_*****.tar`. The `Series Matrix File(s)` usually contains a file called `GSE*****_series_matrix.txt.gz`, where the meta data and processed (usually normalized) data are stored. So in this case you can download and read `GSE*****_series_matrix.txt.gz` file for further analysis. `GSE*****_*****.tar` stores the raw data.

For example, in the project GSE31193, the expression is profiled by Affymetrix Human Genome U133 Plus 2.0 Array. This array is on [GPL570](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL570) (HG-U133_Plus_2) platform. The information for the probes on the microarray is the same if different projects are profiled on the same platform. When you open [GSE31193](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE31193) page and move down to the bottom of the webpage, you will find `Download family` area, where there is a `Series Matrix File(s)`. By clicking [Series Matrix File(s)](ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE31nnn/GSE31193/matrix/) you will be directed to a ftp webpage, where you can download the data by simply click the file name. Unless you really need the raw data in this project, you have already got enough information for most of the gene expression data analysis. If you need the raw data, then simply download the raw data by click `GSE*****_*****.tar` for further analysis.

[GEOquery](https://bioconductor.org/packages/release/bioc/html/GEOquery.html) package gives you an ideal option to deal with this scenario. You don't 
have to download the `Series Matrix File(s)` to get the data.

```r
# install GEOquery package
if (!requireNamespace("GEOquery"))
 BiocManager::install("GEOquery")

library(GEOquery)

# expressionSet object
esets = getGEO(GEO = "GSE31193")
esets

# what you need is the first element in the list
eset = esets[[1]]

# phenodata or metadata
pd = pData(eset)

# expression data
expr = exprs(eset)

# probes information
probes = fData(eset)

# probe information and the meaning of each column
fd = featureData(eset)
```

Of note, GEOquery does not download the raw data automatically, you need to manually download the `GSE*****_*****.tar` for further analysis.

## RNA-seq data
For RNA-seq data, the data are usually stored not exactly the same as for microarray data.

The same thing is that the `Series Matrix File(s)` also contains a file called `GSE*****_series_matrix.txt.gz`. But the difference is that in RNA-seq projects this file stores only the meta data but not the processed/raw expression data. So by using GEOquery package, it does not give you enough information to perform gene expression analysis. The good thing is that GEO database requires the authors of the project to upload either raw or normalized read count data as a supplementary file (usually called `GSE*****_*****.tar`). You still need to manually download the supplementary file and read it into R for further gene expression data analysis.

### Download metadata
For example, [GSE130567](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130567) is an RNA-seq project where the gene expression is profiled on the  [GPL20301](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GPL20301) (Illumina HiSeq 4000) platform. In this project, we need to download series matrix file [`GSE130567_series_matrix.txt.gz`](https://ftp.ncbi.nlm.nih.gov/geo/series/GSE130nnn/GSE130567/matrix/GSE130567_series_matrix.txt.gz) for recovering meta data, and the supplementary file [`GSE130567_RAW.tar`](https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130567&format=file) for recovering gene expression data. Again GEOquery package can deal with the meta data.

```r
# download data using GEOquery package
eset2 = getGEO(GEO = "GSE130567")[[1]]

# phenodata or metadata
pd2 = pData(eset2)
dim(pd2)

# expression data
expr2 = exprs(eset2) # No expression data (0 row)
dim(expr2)

# probes information
probes2 = fData(eset2)
dim(probes2) # No probes data (0 row)

# probe information and the meaning of each column
fd2 = featureData(eset2)
fd2  # No probes data (none)

```

### Download expression data files
Because GEOquery package does not get the expression data, now we need to manually download expression data, which is stored in GSE130567_RAW.tar file (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130567&format=file). First, we download GSE130567_RAW.tar file into a temporary folder, and uncompressed it to obtain the expression file name for each sample (in total 24 samples). 
```r
# temporary folder name and file name
tmpFolder = tempdir()
tmpFile = tempfile(pattern = "GSE130567_", tmpdir = tmpFolder, fileext = ".tar")

# download the raw read count file
url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130567&format=file"
download.file(url, destfile = tmpFile, mode = "wb")

# uncompressed the file into the folder
folder = paste0(tmpFolder, "/GSE130567_count")
untar(tmpFile, exdir = folder)
# get the uncompresed filenames
files = untar(tmpFile, list = TRUE)
# the full paths of the files
filesFull = file.path(folder, files)

# 24 samples
length(filesFull)
filesFull
```

### Read expression files
According to the previous results, there are 24 files compressed in GSE130567_RAW.tar file. We need to one-by-one read the files.

```r
dat = list()
for (file in filesFull){
  accID = gsub(".*(GSM\\d{7}).*", "\\1", file)
  zz = gzfile(file, "rt")
  zzdata = read.csv(zz, header = FALSE, sep = "\t", skip = 4, row.names = 1)
  close(zz)
  zzdata = zzdata[,1, drop = FALSE] # Extract the first numeric column
  colnames(zzdata) = accID
  dat = c(dat, list(zzdata))
}
edata1 = do.call(cbind, dat[1:12])
edata2 = do.call(cbind, dat[13:24])

# The numbers of rows are different
nrow(edata1)
nrow(edata2)
```


The numbers of rows are different, meaning that the pipeline that is used to generate the raw read counts for the first 12 samples are different from the one for the last 12 samples. So we may need to download the fastq files and use the same pipeline to get the same number of rows of raw read counts if we want to compare 24 samples. And this is the objective of this package (SRRDownloader).

### Download fastq files from SRA database

Before downloading fastq files, it's better to know the names of following IDs:

- Platform ID: GPL***** (for example GPL17021)

- GEO Accession ID: GSE***** (for example GSE71165)

- GEO Sample Name ID: GSM******* (for example GSM1828772)

- SRA ? ID: SRS******* (for example SRS1008583)

- SRA Run ID: SRR******* (for example SRR2121685)

- SRA Study ID: SRP****** (for example SRP061381)

- SRA Experiment ID: SRX******* (for example SRX1114423)

- BioSample ID: SAMN******** (for example SAMN03892459)

- BioProject ID: PRJNA****** (for example PRJNA290485)

#### sra-tools

SRRDownloader package uses [`sra-tools`](https://github.com/ncbi/sra-tools/wiki/) to download the fastq files by providing the SRA Run IDs. So the very first step is to install `sra-tools`. Installing `sra-tools` is OS-dependent. 

For Debian and Ubuntu users, you can type: apt install sra-toolkit in your command line to install the toolkit. You can read more about SRA toolkit [here](https://www.ncbi.nlm.nih.gov/books/NBK242621/) and at their [github](https://github.com/ncbi/sra-tools) repo.



Please use the following command to check how to install `sra-tools`.
```

```


#### Entrez Direct

It's optional to use pure command lines to get the work done. It can be achieved 
by using Entrez Direct ([EDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/)) -- E-utilities on the Unix Command Line.


Please use the following R command to check if you have installed or how to install `EDirect`.

```r
installSRRDownloader()
```

If you are using anaconda, to install this package run the following command in a bash terminal:

```console
foo@bar:~$ conda install -c bioconda entrez-direct
```



