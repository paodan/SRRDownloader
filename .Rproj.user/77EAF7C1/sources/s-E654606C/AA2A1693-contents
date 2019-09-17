#' Search SRR accession IDs using projectID
#' @param projectID character, the project accession ID to search in sra database
#' (https://www.ncbi.nlm.nih.gov/sra). This ID can be searched by
#' \code{\link{searchProjectID}} function.
#' @param timeout numeric, the number of seconds to wait before killing the searching process.
#' @return a SRR documents "runinfo" data frame containing 47 columns, which are
#' Run, ReleaseDate, LoadDate, spots, bases, spots_with_mates, avgLength,
#' size_MB, AssemblyName, download_path, Experiment, LibraryName, LibraryStrategy,
#' LibrarySelection, LibrarySource, LibraryLayout, InsertSize, InsertDev,
#' Platform, Model, SRAStudy, BioProject, Study_Pubmed_id, ProjectID, Sample,
#' BioSample, SampleType, TaxID, ScientificName, SampleName, g1k_pop_code, source,
#' g1k_analysis_group, Subject_ID, Sex, Disease, Tumor, Affection_Status, Analyte_Type,
#' Histological_Type, Body_Site, CenterName, Submission, dbgap_study_accession,
#' Consent, RunHash, and ReadHash.\cr
#' Please note RunHash and ReadHash are not MD5SUM of each file.
#' @import funcTools
#' @export
#' @examples {
#' \dontrun{
#' searchSrrID()
#' searchSrrID("PRJNA543132")
#' }
#' }
searchSrrID = function(projectID = "PRJNA540657", timeout = 60){
  format = "runinfo"
  cmd = sprintf('esearch -db sra -query "%s" | efetch -format %s',
                projectID, format)
  system("which esearch")
  res1 = system(cmd, intern = T, timeout = timeout)

  if (length(res1) == 0){
    return(NULL)
  }
  title = strSplit(grep("^Run", res1, value = T), ",")[1,]
  details = as.data.frame(strSplit(grep("^SRR", res1, value = TRUE), ","), stringsAsFactors = F)
  colnames(details) = title

  return(details)
}
