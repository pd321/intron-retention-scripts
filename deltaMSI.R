suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(tidyr))
suppressPackageStartupMessages(library(purrr))

option_list = list(
  make_option(c("-t", "--trt"), action="store", default=NA, type='character',
              help="MSI file for treatment"),
  make_option(c("-c", "--cnt"), action="store", default=NA, type='character',
              help="MSI file for control"),
  make_option(c("-a", "--trtname"), action="store", default="trt", type='character',
              help="Sample name for treatment. Default: [default %default]"),
  make_option(c("-b", "--cntname"), action="store", default="ctrl", type='character',
              help="Sample file for control. Default: [default %default]"),
  make_option(c("-o", "--out"), action="store", default=NA, type='character',
              help="Outfile")
)
opt = parse_args(OptionParser(option_list=option_list))

# Read in the MSI files
trt.df <- readr::read_tsv(opt$trt, guess_max = 100000)
cnt.df <- readr::read_tsv(opt$cnt, guess_max = 100000)

# Merge the trt and cnt files
trt.cnt.merged.df.raw <- dplyr::inner_join(trt.df, cnt.df, by = "ID", suffix = c("_trt", "_ctrl"))

# Apply QC filters
trt.cnt.merged.df <- trt.cnt.merged.df.raw %>% 
  dplyr::filter(A1_trt > 1 | A1_ctrl > 1) %>% # Junctions should have sufficient coverage
  dplyr::filter(A2_trt > 1 | A2_ctrl > 1) %>%
  dplyr::filter(A1_trt + A2_trt > 4 | A1_ctrl + A2_ctrl > 4) %>% # There should be sufficient total coverage at the juntions
  dplyr::filter(CovFrac_trt > 0.95 | CovFrac_ctrl > 0.95) %>% # There should be retention seen across majority of the intron.
  dplyr::mutate(DeltaMSI = MSI_trt - MSI_ctrl) %>% # Calc difference between the msi
  dplyr::mutate(Intron_trt = A1_trt + A2_trt, Intron_ctrl = A1_ctrl + A2_ctrl) # Used for calculating fisher p-val


fisher.df <- trt.cnt.merged.df %>% 
  dplyr::select(c(Exon_trt, Exon_ctrl, Intron_trt, Intron_ctrl))

pvalues = apply(fisher.df,1, function(x){fisher.test( matrix( c(x[1], x[2], x[3], x[4]), nrow = 2 ) )$p.value})

# Add calculated pvals back to df and apply false discovery rate correction
trt.cnt.merged.df = trt.cnt.merged.df %>%
  dplyr::mutate(pval = pvalues) %>%
  dplyr::mutate(padj = p.adjust(p = pval))

# Rename the sample levels to user given ones
colnames(trt.cnt.merged.df) <- colnames(trt.cnt.merged.df) %>% 
  map_chr(~gsub("trt", opt$trtname, .x)) %>% 
  map_chr(~gsub("ctrl", opt$cntname, .x))

# Select only the relevant cols
trt.cnt.merged.df.out <- trt.cnt.merged.df %>% 
  tidyr::separate(col = "ID", into = c("Chrom", "Start", "End", "Strand"), sep ="\\|") %>%
  dplyr::select(Chrom, Start, End, Strand, DeltaMSI, pval, padj,
                dplyr::starts_with("MSI"), dplyr::starts_with("A1"), dplyr::starts_with("A2"), 
                dplyr::starts_with("Exon"), dplyr::starts_with("CovFrac")) %>%
  dplyr::arrange(padj)

# Write out the output
write.table(x=trt.cnt.merged.df.out, file=opt$out, sep = "\t", row.names = FALSE, quote = FALSE)
