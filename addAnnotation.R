suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))

option_list = list(
  make_option(c("-g", "--gene"), action="store", default=NA, type='character',
              help="intron bed file with genename for each intron e.g. introns.bed"),
  make_option(c("-t", "--type"), action="store", default=NA, type='character',
              help="File with the intron types i.e. U2/U12 e.g. intron_type.xls"),
  make_option(c("-i", "--geneinfo"), action="store", default=NA, type='character',
              help="File with the information for the genes e.g. geneinfo.tsv"),
  make_option(c("-d", "--deltamsi"), action="store", default=NA, type='character',
              help="File giving delta msi values"),
  make_option(c("-o", "--out"), action="store", default=NA, type='character',
              help="File to store the merged annotations")
)
opt = parse_args(OptionParser(option_list=option_list))

processGenes <- function(geneFile, geneInfoFile){
  
  gene.df <- readr::read_tsv(geneFile, col_names = c("Chrom", "Start", "End", "GeneNames", "Score", "Strand"))
  
  gene.df.mutated <- gene.df %>% 
    dplyr::mutate(ID = paste(Chrom, Start, End, Strand, sep = "|")) %>% # Adding ID for join
    dplyr::select(ID, GeneNames) %>% 
    tidyr::separate_rows(GeneNames, sep = ",") %>% # Split by comma so indivudial genes can be mapped to names
    dplyr::mutate(GeneNames = gsub("\\..*", "", GeneNames)) # Remove version numbers makes the mapping easier

  gene.to.name.map <- readr::read_tsv(geneInfoFile) %>%
    dplyr::mutate(Ens_id = gsub("\\..*", "", Ens_id)) %>% # Remove version numbers makes the mapping easier
    dplyr::select(Ens_id, GeneName)
  
  gene.df.name <- gene.df.mutated %>%
    dplyr::left_join(gene.to.name.map, by = c("GeneNames" = "Ens_id")) %>% # map id to gene name
    dplyr::select(ID, GeneName) %>%
    dplyr::group_by(ID) %>%
    dplyr::summarise(Gene = paste(GeneName, collapse = ",")) # Collapse back to the way it was earlier
    
  return(gene.df.name)
}

processType <- function(typeFile){
  type.df <- readr::read_tsv(typeFile)
  
  type.df.mutate <- type.df %>%
    dplyr::mutate(ID = paste(Chrom, Start, End, Strand, sep = "|")) %>%
    dplyr::select(ID, Type, SubType, Confidence)
  
  return (type.df.mutate)
}


deltaMSI.df <- readr::read_tsv(opt$deltamsi)
deltaMSI.df.mutate <- deltaMSI.df %>%
  dplyr::mutate(ID = paste(Chrom, Start, End, Strand, sep = "|"))

out.df <- Reduce(function(x, y) merge(x, y, by = "ID"), list(deltaMSI.df.mutate, processGenes(opt$gene, opt$geneinfo), processType(opt$type)))
out.df <- out.df %>%
  dplyr::select(-ID) %>%
  dplyr::select(Chrom, Start, End, Strand, DeltaMSI, pval, padj , Gene, Type, SubType, Confidence, everything()) %>%
  dplyr::mutate(Start = Start + 1) %>% # Adjusting the start value as we have been using the bed based start values
  dplyr::arrange(padj)

write.table(x=out.df, file=opt$out, sep = "\t", row.names = FALSE, quote = FALSE)
