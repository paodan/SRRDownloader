#' Download fastq files from SRA database by providing SRR accession IDs.
#' @param srrIDs character, the SRR accession ID to download. This ID can be searched by
#' \code{\link{searchSrrID}} function.
#' @param timeout numeric, the number of seconds to wait before killing the
#' downloading process.
#' @param OutDir character, the path to save downloaded fastq files.
#' The default is current directory ("./").
#' @param multipleDownload integer, the number of downloading jobs to run in each batch.
#' After each batch, the program will wait for 100 seconds to initiate next batch.
#' multipleDownload must be 1 if the Rstudio version < 1.2.
#' @return Downloaded reads if multipleDownload = 1, or the temporary job file names if
#' multipleDownload > 1.
#'
#' @details
#' Frequently Used Options:\cr
#' General:\cr
#' -h 	| 	--help 	Displays ALL options, general usage, and version information.\cr
#' -V 	| 	--version 	Display the version of the program.\cr
#' \cr
#' Data formatting:\cr
#' --split-files 	Dump each read into separate file. Files will receive suffix corresponding to read number.\cr
#' --split-spot 	Split spots into individual reads.\cr
#' --fasta <[line width]> 	FASTA only, no qualities. Optional line wrap width (set to zero for no wrapping).\cr
#' -I 	| 	--readids 	Append read id after spot id as 'accession.spot.readid' on defline.\cr
#' -F 	| 	--origfmt 	Defline contains only original sequence name.\cr
#' -C 	| 	--dumpcs <[cskey]> 	Formats sequence using color space (default for SOLiD). "cskey" may be specified for translation.\cr
#' -B 	| 	--dumpbase 	Formats sequence using base space (default for other than SOLiD).\cr
#' -Q 	| 	--offset <integer> 	Offset to use for ASCII quality scores. Default is 33 ("!").\cr
#' \cr
#' Filtering:\cr
#' -N 	| 	--minSpotId <rowid> 	Minimum spot id to be dumped. Use with "X" to dump a range.\cr
#' -X 	| 	--maxSpotId <rowid> 	Maximum spot id to be dumped. Use with "N" to dump a range.\cr
#' -M 	| 	--minReadLen <len> 	Filter by sequence length >= <len>\cr
#' --skip-technical 	Dump only biological reads.\cr
#' --aligned 	Dump only aligned sequences. Aligned datasets only; see sra-stat.\cr
#' --unaligned 	Dump only unaligned sequences. Will dump all for unaligned datasets.\cr
#' \cr
#' Workflow and piping:\cr
#' -O 	| 	--outdir <path> 	Output directory, default is current working directory ('.').\cr
#' -Z 	| 	--stdout 	Output to stdout, all split data become joined into single stream.\cr
#' --gzip 	Compress output using gzip.\cr
#' --bzip2 	Compress output using bzip2.\cr
#'
#' @import funcTools
#' @import rstudioapi
#' @export
#' @examples {
#' \dontrun{
#' # Download one by one
#' downloadSrr(srrIDs = c("SRR9063863", "SRR9063864"), OutDir = "./down")
#' # Download 2 files in one batch and wait for 100 seconds and then download the second batch ...
#' downloadSrr(srrIDs = c("SRR9063863", "SRR9063864"), OutDir = "./down", multipleDownload = 2)
#'
#' # Download control samples and interferon treated samples in PRJNA540657 project.
#' x = searchSrrID("PRJNA540657")
#' x = subset(x, SampleName %in% c("GSM3743639", "GSM3743640", "GSM3743641",
#'            "GSM3743645", "GSM3743646", "GSM3743647"))$Run
#' downloadSrr(srrIDs = x, OutDir = "./down", multipleDownload = 3)
#' }
#' }
downloadSrr = function(srrIDs = c("SRR9063863", "SRR9063864"), timeout = 3600*12,
                       OutDir = "./", multipleDownload = 1){

  stopifnot(length(srrIDs) > 0)
  stopifnot((timeout<-as.integer(timeout)) > 0)
  stopifnot((multipleDownload<-as.integer(multipleDownload)) > 0)
  if (multipleDownload > 1){
    if (!"jobRunScript" %in% ls("package:rstudioapi")){
      multipleDownload = 1
    }
  }
  if(!dir.exists(OutDir)){
    dir.create(OutDir, recursive = TRUE)
  }
  if (multipleDownload == 1 | length(srrIDs) == 1){
    cmd = sprintf('fastq-dump %s --split-files -O "%s" --gzip',
                  paste0(srrIDs, collapse = " "), normalizePath(OutDir))

    system("which fastq-dump")
    system.time({
      res1 = system(cmd, intern = T, timeout = timeout)
    })
  } else {
    nBatch = ceiling(length(srrIDs) / multipleDownload)
    tmpFiles = c()
    for (ni in 1:nBatch){
      m = ((ni -1)*multipleDownload+1) : min(length(srrIDs), (ni)*multipleDownload)
      for (mi in m){
        tmpFile = tempfile(pattern = srrIDs[mi], fileext = ".R", tmpdir = tempdir())
        script = sprintf("echo 'sink(\"%s\"); system(\"fastq-dump %s --split-files -O %s --gzip\", intern = FALSE, timeout = %s)'",
                         paste0(tmpFile, ".log"), srrIDs[mi], normalizePath(OutDir), timeout)
        system(paste0(script, " > ", tmpFile))
        rstudioapi::jobRunScript(path = tmpFile, importEnv = TRUE)
        tmpFiles = c(tmpFiles, tmpFile)
      }
      Sys.sleep(100)
    }
    res1 = tmpFiles
  }
  return(res1)
}
