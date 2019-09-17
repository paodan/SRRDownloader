#' Search SRR project accession ID using query
#' @param query character, the title or keywords to search in BioProject database
#' (https://www.ncbi.nlm.nih.gov/bioproject)
#' @param timeout numeric, the number of seconds to wait before killing the searching process.
#'
#' @details
#' Format Selection\cr
#'\cr
#' -format        Format of record or report\cr
#' -mode          text, xml, asn.1, json\cr
#' -style         withparts, conwithfeat\cr
#'\cr
#' Direct Record Selection\cr
#'\cr
#' -db            Database name\cr
#' -id            Unique identifier or accession number\cr
#'\cr
#' Sequence Range\cr
#'\cr
#' -seq_start     First sequence position to retrieve\cr
#' -seq_stop      Last sequence position to retrieve\cr
#' -strand        Strand of DNA to retrieve\cr
#'\cr
#' Gene Range\cr
#'\cr
#' -chr_start     Sequence range from 0-based coordinates\cr
#' -chr_stop        in gene docsum GenomicInfoType object\cr
#'\cr
#' Sequence Flags\cr
#'\cr
#' -complexity    0 = default, 1 = bioseq, 3 = nuc-prot set\cr
#' -extend        Extend sequence retrieval in both directions\cr
#' -extrafeat     Bit flag specifying extra features\cr
#'\cr
#' Miscellaneous\cr
#'\cr
#' -raw           Skip database-specific XML modifications\cr
#' -json          Convert adjusted XML output to JSON\cr
#'\cr
#' Format Examples\cr
#'\cr
#' -db            -format            -mode    Report Type\cr
#' ___            _______            _____    ___________\cr
#'\cr
#' (all)\cr
#' docsum                      DocumentSummarySet XML\cr
#' docsum             json     DocumentSummarySet JSON\cr
#' full                        Same as native except for mesh\cr
#' uid                         Unique Identifier List\cr
#' url                         Entrez URL\cr
#' xml                         Same as -format full -mode xml\cr
#'\cr
#' bioproject\cr
#' native                      BioProject Report\cr
#' native             xml      RecordSet XML\cr
#'\cr
#' biosample\cr
#' native                      BioSample Report\cr
#' native             xml      BioSampleSet XML\cr
#'\cr
#' biosystems\cr
#' native             xml      Sys-set XML\cr
#'\cr
#' gds\cr
#' native             xml      RecordSet XML\cr
#' summary                     Summary\cr
#'\cr
#' gene\cr
#' full_report                 Detailed Report\cr
#' gene_table                  Gene Table\cr
#' native                      Gene Report\cr
#' native             asn.1    Entrezgene ASN.1\cr
#' native             xml      Entrezgene-Set XML\cr
#' tabular                     Tabular Report\cr
#'\cr
#' homologene\cr
#' alignmentscores             Alignment Scores\cr
#' fasta                       FASTA\cr
#' homologene                  Homologene Report\cr
#' native                      Homologene List\cr
#' native             asn.1    HG-Entry ASN.1\cr
#' native             xml      Entrez-Homologene-Set XML\cr
#'\cr
#' mesh\cr
#' full                        Full Record\cr
#' native                      MeSH Report\cr
#' native             xml      RecordSet XML\cr
#'\cr
#' nlmcatalog\cr
#' native                      Full Record\cr
#' native             xml      NLMCatalogRecordSet XML\cr
#'\cr
#' pmc\cr
#' medline                     MEDLINE\cr
#' native             xml      pmc-articleset XML\cr
#'\cr
#' pubmed\cr
#' abstract                    Abstract\cr
#' medline                     MEDLINE\cr
#' native             asn.1    Pubmed-entry ASN.1\cr
#' native             xml      PubmedArticleSet XML\cr
#'\cr
#' (sequences)\cr
#' acc                         Accession Number\cr
#' est                         EST Report\cr
#' fasta                       FASTA\cr
#' fasta              xml      TinySeq XML\cr
#' fasta_cds_aa                FASTA of CDS Products\cr
#' fasta_cds_na                FASTA of Coding Regions\cr
#' ft                          Feature Table\cr
#' gb                          GenBank Flatfile\cr
#' gb                 xml      GBSet XML\cr
#' gbc                xml      INSDSet XML\cr
#' gene_fasta                  FASTA of Gene\cr
#' gp                          GenPept Flatfile\cr
#' gp                 xml      GBSet XML\cr
#' gpc                xml      INSDSet XML\cr
#' gss                         GSS Report\cr
#' ipg                         Identical Protein Report\cr
#' ipg                xml      IPGReportSet XML\cr
#' native             text     Seq-entry ASN.1\cr
#' native             xml      Bioseq-set XML\cr
#' seqid                       Seq-id ASN.1\cr
#'\cr
#' snp\cr
#' chr                         Chromosome Report\cr
#' docset                      Summary\cr
#' fasta                       FASTA\cr
#' flt                         Flat File\cr
#' native             asn.1    Rs ASN.1\cr
#' native             xml      ExchangeSet XML\cr
#' rsr                         RS Cluster Report\cr
#' ssexemplar                  SS Exemplar List\cr
#'\cr
#' sra\cr
#' native             xml      EXPERIMENT_PACKAGE_SET XML\cr
#' runinfo            xml      SraRunInfo XML\cr
#'\cr
#' structure\cr
#' mmdb                        Ncbi-mime-asn1 strucseq ASN.1\cr
#' native                      MMDB Report\cr
#' native             xml      RecordSet XML\cr
#'\cr
#' taxonomy\cr
#' native                      Taxonomy List\cr
#' native             xml      TaxaSet XML\cr
#' @export
#' @examples {
#' \dontrun{
#' searchProjectID()
#' searchProjectID("IFN-γ T cells")
#' }
#' }
searchProjectID = function(query = "IFN-γ selectively suppresses a subset of TLR4-activated genes and enhancers to potentiate macrophage activation",
                           timeout = 60){
  format = "native"
  cmd = sprintf('esearch -db bioproject -query "%s" | efetch -format %s',
                query, format)
  system("which esearch")
  res1 = system(cmd, intern = TRUE, timeout = timeout)
  index = which(res1 == "")
  len = length(index)/2
  if(len == 0){
    return(NULL)
  }
  res = as.data.frame(matrix("", nrow = len, ncol = 6,
                             dimnames = list(1:len, c("index", "title", "organism", "accession", "ID", "note"))),
                      stringsAsFactors = FALSE)
  .fillBlank = function(x, start, end){
    if (length(x) == 0){
      return("")
    } else {
      return(substr(x, start, end))
    }
  }
  for (ni in 1 : len){
    mi = seq(1, len*2, by = 2)[ni]
    tmp = res1[index[mi]: index[mi + 1]]
    res[ni, "index"] = ni
    organism = grep("Organism: ", tmp, value = TRUE)
    accession = grep("BioProject Accession: ", tmp, value = TRUE)
    ID = grep("ID: ", tmp, value = TRUE)
    title = sub("^[0-9]+\\. ", "", tmp[2])

    res[ni, "organism"] = .fillBlank(organism, 11, 1000)
    res[ni, "accession"] = .fillBlank(accession, 23, 100)
    res[ni, "ID"] = .fillBlank(ID, 5, 100)
    res[ni, "note"] = tmp[length(tmp)]
    res[ni, "title"] = title
  }
  return(res)
}
