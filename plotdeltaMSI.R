suppressPackageStartupMessages(library(optparse))
suppressPackageStartupMessages(library(dplyr))
suppressPackageStartupMessages(library(readr))
suppressPackageStartupMessages(library(purrr))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(cowplot))

option_list = list(
  make_option(c("-i", "--deltaMSI"), action="store", default=NA, type='character',
              help="deltaMSI output file from deltaMSI.R"),
  make_option(c("-o", "--outdot"), action="store", default=NA, type='character',
              help="Outfile for dot plot"),
  make_option(c("-n", "--outdensity"), action="store", default=NA, type='character',
              help="Outfile for density plot")
)
opt = parse_args(OptionParser(option_list=option_list))

# Read in the deltaMSI file
delta.msi.df <- readr::read_tsv(opt$deltaMSI, guess_max = 100000)

# Add negLog10Pval
delta.msi.df <- delta.msi.df %>% 
  dplyr::mutate(negLog10Pval = -log10(pval))

# Set int_type as factor
delta.msi.df$int_type <- factor(delta.msi.df$Type, levels = c("U12", "U2"))

# Remove outliers
pval_percentile <- stats::quantile(delta.msi.df$negLog10Pval, probs = c(0.9999))
delta_msi_percentile <- stats::quantile(delta.msi.df$DeltaMSI, probs = c(0.001,0.999))

delta.msi.df <- delta.msi.df %>% 
  dplyr::filter(negLog10Pval < pval_percentile) %>% 
  dplyr::filter(DeltaMSI > delta_msi_percentile[1], DeltaMSI < delta_msi_percentile[2])

# Plot Dot Plot
dot_plot <- ggplot2::ggplot(data = delta.msi.df, 
                            mapping = aes(x = DeltaMSI, y = negLog10Pval, label = Gene, color = Type)) +
  scale_colour_manual(name = "IntronType", values = c("red", "grey")) + 
  geom_point(size = 1.5) + 
  geom_point(data = delta.msi.df %>% 
    dplyr::filter(int_type == "U12"), color = "red",  alpha = 0.9, size = 2) +
  ylab("-log10(pval)") + 
  xlab(expression(paste(Delta,"MSI")))+
  geom_vline(xintercept = 0, linetype = "dashed") + 
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  theme(axis.text=element_text(size=20), axis.title.x =element_text(size=20), 
        axis.title.y=element_text(size=20), legend.position="none")


# Save Dot Plot
cowplot::save_plot(filename = opt$outdot, plot = dot_plot, base_height = 10, base_width = 10)

# Plot Density Plot
density_plot <- ggplot(delta.msi.df, aes(colour = int_type)) + 
  stat_density(aes(DeltaMSI), geom="line",position="identity") + 
  scale_colour_manual(name = "IntronType", values = c("red", "black"))  + 
  geom_density(data = delta.msi.df %>% 
    dplyr::filter(int_type == "U2"), aes(DeltaMSI), color = "black") + 
  geom_density(data = delta.msi.df %>% 
    dplyr::filter(int_type == "U12"), aes(DeltaMSI), color = "red")

# Save Density plot
cowplot::save_plot(filename = opt$outdensity, plot = density_plot, base_height = 10, base_width = 10)



