# Download gene expression data from GEO database
Author: Weiyang Tao

Email: weiyangtao1513@gmail.com

Date: 2023-08-02

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

### Download expression data files and read expression files
Because GEOquery package does not get the expression data, now we need to manually download expression data, which is stored in GSE130567_RAW.tar file (https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE130567&format=file). First, we download GSE130567_RAW.tar file into a temporary folder, and uncompressed it to obtain the expression file name for each sample (in total 24 samples). Then we need to one-by-one read the files. The `SRRDownloader` package provide a function `getSupplementary` to download and read the supplementary file in one go.

```r
library(SRRDownloader)
supp = getSupplementaryData("GSE130567", header = FALSE)
dat = supp$data
names(dat) = sapply(strsplit(names(dat), "_"), "[", 1)

edata1 = do.call(cbind, dat[1:12])
edata2 = do.call(cbind, dat[13:24])

# The numbers of rows are different
nrow(edata1)
nrow(edata2)
```

Alternatively, you can read it manually, but it is much more complicated.
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

# Read the expression data
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


The numbers of rows of `edata1` and `edata2` are different, meaning that the pipeline that is used to generate the raw read counts in the original paper for the first 12 samples are different from the one for the last 12 samples. So we may need to download the fastq files and use the same pipeline to get the same number of rows of raw read counts if we want to compare 24 samples. And this is the objective of this package (SRRDownloader).

### Download fastq files from SRA database

Before downloading fastq files, it's better to know the names of following IDs:

- Platform ID: `GPL*****` (for example GPL17021)

- GEO Accession ID: `GSE*****` (for example GSE71165)

- GEO Sample Name ID: `GSM*******` (for example GSM1828772)

- SRA ? ID: `SRS*******` (for example SRS1008583)

- SRA Run ID: `SRR*******` (for example SRR2121685)

- SRA Study ID: `SRP******` (for example SRP061381)

- SRA Experiment ID: `SRX*******` (for example SRX1114423)

- BioSample ID: `SAMN********` (for example SAMN03892459)

- BioProject ID: `PRJNA******` (for example PRJNA290485)


#### How to install sra-tools

`SRRDownloader` package uses [`sra-tools`](https://github.com/ncbi/sra-tools/wiki/) to download the fastq files by providing the SRA Run IDs (`SRR*******`, this is where the name of the package comes from). So the very first step is to install `sra-tools`. Installing `sra-tools` is OS-dependent. You can read more about SRA toolkit [here](https://www.ncbi.nlm.nih.gov/books/NBK242621/) and at their [github](https://github.com/ncbi/sra-tools) repo.

You can follow these steps to install the toolkit. 

##### - For Windows

1. Download the sratoolkit.current-win64.zip file from this [link](https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-win64.zip)

2. Extract it to your desktop, for example

3. Open a command shell, for example Start/Run `cmd.exe`

4. `cd` to the directory you extracted the zip file to, for example Desktop

5. `cd bin`

6. Proceed to the [Quick Configuration Guide](https://github.com/ncbi/sra-tools/wiki/03.-Quick-Toolkit-Configuration)

7. Test that the toolkit is functional:

```console
fastq-dump --stdout -X 2 SRR390728
```

Within a few seconds, the command should produce this exact output (and nothing else):

```console
Read 2 spots for SRR390728
Written 2 spots for SRR390728
@SRR390728.1 1 length=72
CATTCTTCACGTAGTTCTCGAGCCTTGGTTTTCAGCGATGGAGAATGACTTTGACAAGCTGAGAGAAGNTNC
+SRR390728.1 1 length=72
;;;;;;;;;;;;;;;;;;;;;;;;;;;9;;665142;;;;;;;;;;;;;;;;;;;;;;;;;;;;;96&&&&(
@SRR390728.2 2 length=72
AAGTAGGTCTCGTCTGTGTTTTCTACGAGCTTGTGTTCCAGCTGACCCACTCCCTGGGTGGGGGGACTGGGT
+SRR390728.2 2 length=72
;;;;;;;;;;;;;;;;;4;;;;3;393.1+4&&5&&;;;;;;;;;;;;;;;;;;;;;<9;<;;;;;464262
```


##### - For Debian and Ubuntu

Use the following command:

```console
apt install sra-toolkit
```

##### - For MAC OS X

Use the following command:

```console
# 1. Fetch the tar file from the canonical location at NCBI:
curl --output sratoolkit.tar.gz https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/current/sratoolkit.current-mac64.tar.gz
# 2. Extract the contents of the tar file:
tar -vxzf sratoolkit.tar.gz
# 3. For convenience (and to show you where the binaries are) append the path to the binaries to your PATH environment variable (the version number may be different), and add it to ~/.bash_profile file:
export PATH=$PATH:$PWD/sratoolkit.3.0.0-mac64/bin
# 4. Verify that the binaries will be found by the shell:
which fastq-dump
```
Follow the step 6-7 in **For Windows** section.


<br>

#### How to install Entrez Direct tool

Getting SRA Run IDs requires you searching in GEO database. You might ask if there is any tool that can search the SRA Run IDs using command lines. The answer is YES, you can search SRA Run IDs
by using Entrez Direct ([EDirect](https://www.ncbi.nlm.nih.gov/books/NBK179288/)) -- E-utilities on the Unix Command Line. But this is optional for the use of `SRRDownloader` package.


Please use the following R command to check if you have installed or how to install `EDirect`.

```r
installSRRDownloader()
```


EDirect will run on Unix and Macintosh computers, and under the Cygwin Unix-emulation environment on Windows PCs. To install the EDirect software, open a terminal window and execute one of the following two commands:
```console
sh -c "$(curl -fsSL https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh)"

sh -c "$(wget -q https://ftp.ncbi.nlm.nih.gov/entrez/entrezdirect/install-edirect.sh -O -)"
```
This will download a number of scripts and several precompiled programs into an "edirect" folder in the user's home directory. It may then print an additional command for updating the PATH environment variable in the user's configuration file. The editing instructions will look something like:

```
echo "export PATH=\$HOME/edirect:\$PATH" >> $HOME/.bash_profile
```

As a convenience, the installation process ends by offering to run the PATH update command for you. Answer "y" and press the Return key if you want it run. If the PATH is already set correctly, or if you prefer to make any editing changes manually, just press Return.

One installation is complete, run:

```console
export PATH=${HOME}/edirect:${PATH}
```
to set the PATH for the current terminal session.

<br>
*Or if you are using anaconda, to install this package by simply running the following command in a bash terminal:*

```console
conda install -c bioconda entrez-direct
```

<br>
Now you have finished installing tools for downloading and searching. Let's download an example fastq file into the current folder using SRA Run ID [`SRR8996103`](https://trace.ncbi.nlm.nih.gov/Traces/index.html?view=run_browser&acc=SRR8996103&display=metadata) for GEO Sample Name `GSM3743639` in GEO Accession [`GSE130567`](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE130567).

```r
downloadSrr(srrIDs = c("SRR9063863"), timeout = 3600 *12, OutDir = "./", multipleDownload = 1)
```


