#' Submit data files to GEO database by sftp in a fast way, using ncft.
#' @description The detailed submission process can be found in GEO File Transfer Protocol:
#' https://www.ncbi.nlm.nih.gov/geo/info/submissionftp.html.
#' @param fileNames a vector of filenames with full path name to submit
#' @param destDir the directory in GEO FTP server. This directory must be created in advance
#' (in shell):
#' sftp geo@sftp-private.ncbi.nlm.nih.gov \cr
#' password: 33%9uyj_fCh?M16H \cr
#' mkdir destDir \cr
#' \cr
#' More detailes can be found in GEO File Transfer Protocol\cr
#' @param geoFTP GEO FTP server address. Default is "ftp-private.ncbi.nlm.nih.gov".
#' This password can be found in GEO File Transfer Protocol.
#' @param geoPassword GEO FTP server password. Default is "33%9uyj_fCh?M16H".
#' This password can be found in GEO File Transfer Protocol.
#' @param nPerBatch the number of uploading process in each batch. Default is 10.
#' @param sleepSecond the time (in second) to wait after each batch. Default is 100 (seconds).
#' @export
#' @examples {
#' \dontrun{
#' # Upload files in "./tmp" folder to "ftpDir"
#' files = dir("./tmp", full.names = T)
#' fastGeoSubmission(files, "ftpDir")
#' }
#' }
fastGeoSubmission = function(fileNames, destDir, geoFTP = "ftp-private.ncbi.nlm.nih.gov",
                             geoPassword = "33%9uyj_fCh?M16H", nPerBatch = 10, sleepSecond = 100){
  len = length(fileNames)
  if(is.null(len) || len<1){
    stop("fileNames must be not empty.")
  }

  if(length(system("which ncftpput", intern = T)) < 1){
    stop("Please install ncft by following command in bash:\n\n (sudo) apt-get install ncftp\n\n")
  }

  if (is.null(destDir) || length(destDir)!=1){
    stop("destDir must be not empty and of length one.")
  }

  for (ni in seq(1, len, by = nPerBatch)){
    for (mi in fileNames[ni:min(len, ni+nPerBatch-1)]){
      rfile = paste0(tempdir(), "/", basename(mi), ".R")
      cmd = paste0('system("ncftpput -F -R -z -u geo -p \\"', geoPassword, '\\" ',
                   geoFTP, ' ', destDir, ' ', mi, '")')
      cat(cmd, "\n")

      fileConn<-file(rfile)
      writeLines(cmd, fileConn)
      close(fileConn)

      rstudioapi::jobRunScript(rfile)
    }
    Sys.sleep(sleepSecond)
  }
  return(invisible(NULL))
}


