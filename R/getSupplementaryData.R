#' This is to read the supplementary .tar file in GEO database. In .tar file, there are
#' usually several .txt.gz (.csv.gz, .tab.gz) files or one .txt.gz (.csv.gz, .tab.gz) file.
#' @param tarfile the name of supplementary .tar file
#' @param sep the field separator character. Values on each line of the file are
#' separated by this character. The default is `\\t`.
#' @param header a logical value indicating whether the file contains the names
#' of the variables as its first line. If missing, the value is determined from
#' the file format: header is set to TRUE if and only if the first row contains
#' one fewer field than the number of columns. The default is FALSE.
#' @return a list of data frame coming from the .gz files in the .tar file.
#' @export
readTarFile = function(tarfile, sep = "\t", header = FALSE){
  # temporary folder name and file name
  tmpFolder = tempdir()
  folder = tempfile(pattern = "Out_", tmpdir = tmpFolder)

  # uncompressed the file into the folder
  untar(tarfile, exdir = folder)
  # get the uncompresed filenames
  files = untar(tarfile, list = TRUE)
  # the full paths of the files
  filesFull = file.path(folder, files)

  res = lapply(filesFull, function(tmpFile){
    print(tmpFile)
    if(grep(".txt.gz$|.csv.gz$|.tab.gz$", tmpFile)){
      x = read.csv(gzfile(tmpFile), sep = sep, header = header)
    } else {
      stop("Unknown .gz file type in .tar file")
    }
    x
  })
  names(res) = basename(filesFull)
  return(res)
}

#' Download Supplementary file and get the data using GEO accession ID
#' @param accession GEO accession ID, the default is "GSE121213", sep = `\\t`, header = TRUE, suppFile = NULL
#' @param sep the field separator character. Values on each line of the file are
#' separated by this character. The default is `\\t`.
#' @param header a logical value indicating whether the file contains the names
#' of the variables as its first line. If missing, the value is determined from
#' the file format: header is set to TRUE if and only if the first row contains
#' one fewer field than the number of columns. The default is FALSE.
#' @param suppFile supplementary file name. In some of the GEO projects, there are
#' multiple supplementary files provided. In this case, `getSupplementaryData`
#' function by default reads the first file and shows a warning. If a supplementary
#' file name is specified by this parameter, `getSupplementaryData` will read the
#' data using this file name without showing a warning.
#' is specified,
#' @return a list, comprised of `data`, file `extension`, `accession`, `filename`,
#' and `url`.
#' @import rvest
#' @importFrom tools file_ext
#' @import xml2
#' @export
#' @examples
#' \dontrun{
#' # Get the suplementary data for accession ID: GSE121213
#' x = getSupplementaryData()
#' x
#' x$data
#'
#' # Get the suplementary data for accession ID: GSE121212
#' x = getSupplementaryData("GSE121212")
#' x
#' x$data
#'
#' # Get the suplementary data for accession ID: GSE130567
#' # This suplementary file is a .tar file, which contains 24 .txt.gz files.
#' x = getSupplementaryData("GSE130567", header = FALSE)
#' x
#' x$data
#'
#' # You will get a warning because of multiple supplementary files
#' x = getSupplementaryData("GSE69626")
#' x$data
#' x$filename
#'
#' # Specify a file name to avoid the warning
#' x = getSupplementaryData("GSE69626", suppFile = "GSE69626_RPKMGEO.txt.gz")
#' x$data
#' x$filename
#' }
getSupplementaryData = function(accession = "GSE121213", sep = "\t", header = TRUE, suppFile = NULL){
  cat(accession, ": ")
  url = paste0("https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=", accession)
  html <- read_html(url)
  nodes = html_nodes(html, xpath="//a")
  id = which(html_text(nodes) == "(http)")

  if(length(id) == 0) {
    stop("Supplementary data files not provided\n")
  } else if(length(id) > 1){
    # table text
    tblNode = head(tail(html_nodes(html, "table"), 2), 1)
    tblText = html_table(tblNode, fill = T)[[1]][-1, ]
    whichRow = grep("http", tblText[,3])
    filenamesAll = tblText[whichRow,1]
    # file urls
    html_attrs(html_nodes(tblNode, xpath = "//tr[(((count(preceding-sibling::*) + 1) = 2) and parent::*)]//td[(((count(preceding-sibling::*) + 1) = 1) and parent::*)]"))
    suppUrls0 = html_nodes(tblNode, "a") %>% xml_attr("href")
    suppUrls = grep("/geo/download", suppUrls0, value = TRUE)
    suppUrls = sapply(suppUrls, URLdecode)

    if(!is.null(suppFile)){
      stopifnot(length(suppFile) == 1)
      idx = which(filenamesAll == suppFile)
      len = length(idx)
      if(len == 0){
        stop("No this file: ", suppFile)
      } else if(len == 1){
        fileURL = suppUrls[idx]
        filename = filenamesAll[idx]
      } else {
        stop("multiple files are found")
      }
    } else {
      fileURL = suppUrls[1]
      filename = filenamesAll[1]
      warning(length(suppUrls), " files (", paste(filenamesAll, collapse = ", "),
              ") are found, but only the first file will be read: ", filename,
              "\nIf you want to read another file, please specify it by `suppFile` parameter.")
    }
  } else {
    fileURL = html_attrs(nodes[id])[[1]]
    nodes = html_nodes(html, xpath="//tr[(((count(preceding-sibling::*) + 1) = 2) and parent::*)]//td[(((count(preceding-sibling::*) + 1) = 1) and parent::*)]")
    filename = html_text(tail(nodes, 1))
    cat(filename, "\n")
  }

  fileURL = paste0("https://www.ncbi.nlm.nih.gov", fileURL)

  # temporary folder name and file name
  tmpFolder = tempdir()
  tmpFile = paste0(tmpFolder, "/", filename)
  cat("Downloading", tmpFile, "\n")
  s = download.file(fileURL, destfile = tmpFile, mode = "wb")
  if(s != 0) stop("Failed to download ", tmpFile)
  ext = file_ext(tmpFile)

  cat("Reading", tmpFile, "\n")
  if(ext == "gz"){
    if(grep(".txt.gz$|.csv.gz$|.tab.gz$", tmpFile)){
      x = read.csv(gzfile(tmpFile), sep = sep, header = header)
    } else {
      stop("Unknown supplementary .gz file type")
    }
  } else if(ext == "tar"){
    x = readTarFile(tmpFile, sep = sep, header = header)
  }
  cat("Done\n")
  return(list(data = x, extension = ext, accession = accession,
              filename = filename, url = URLdecode(fileURL)))
}


