#' Instruction to install prerequisite softwares.
#' @export
#' @examples {
#' \donotrun{
#' installSRRDownloader()
#' }
#' }
installSRRDownloader = function(){
  cmd1 = system("which esearch", intern = TRUE)
  cmd2 = system("which fastq-dump", intern = TRUE)

  res = c()
  if(length(cmd1) == 0){
    cat("Open a Terminal to install ncbi-entrez-direct by:\n")
    message("sudo apt-get install ncbi-entrez-direct\n")
    res = c(res, "sudo apt-get install ncbi-entrez-direct")
  }

  if(length(cmd2) == 0){
    cat("Open a Terminal to install sra-toolkit by:\n")
    message("sudo apt install sra-toolkit\n")
    res = c(res, "sudo apt install sra-toolkit")
  }

  if(length(cmd1) != 0 && length(cmd2) != 0){
    cat("Prerequisite softwares have been installed.\n")
  }
  return(invisible(res))
}
