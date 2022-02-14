library(ggplot2)
library(tidyverse)
library(reshape2)
library(optparse)

option_list <- list(
    make_option(c("-i", "--input"), type = "character",
        help = "Input table. Multiple tables sebarated by comma.",
        dest = "input"),
    make_option(c("-x", "--xname"), type = "character", default = "Score",
        help = "X-lab name."),
    make_option(c("-y", "--yname"), type = "character", default = "Value",
        help = "Y-lab dname."),
    make_option(c("-o", "--output"), type = "character",
        help = "Output file of the plot")
)

arguments <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE) # nolint

opt <- arguments$options

if (!"input" %in% names(opt)) {
  stop("--input parameter must be provided. See script usage (--help)")
}
if (!"output" %in% names(opt)) {
  stop("--output parameter must be provided. See script usage (--help)")
}



colours_nice <- c("#1b9e77","#d95f02","#7570b3","#e6ab02","#e7298a","#a6761d","#000000") # nolint
size_text <- 25;
size_legend <- 25;
size_line <- 1;
size_geom_line <- 2;
standard_style <- theme_bw() + theme(axis.text = element_text(colour = "black",size=size_text), axis.title = element_text(colour = "black",size=size_text), plot.title = element_text(size = size_text, face="bold"), panel.grid.major = element_blank() , panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks = element_line(colour = "black", size=1), axis.line = element_line(colour = "black", size=size_line), legend.key =  element_blank(), legend.text = element_text(size=size_legend), legend.position="top", legend.direction="horizontal", legend.key.size = unit(2, 'lines'), legend.title=element_blank())+ theme(axis.line.x = element_line(color="black", size = size_line), axis.line.y = element_line(color="black", size = size_line)) # nolint


metrics_per_threshold_tsv <- read_delim(opt$input, delim = "\t",
    escape_double = FALSE, trim_ws = TRUE)
maxf1 <- metrics_per_threshold_tsv %>% select(Threshold, `F1-score`) %>% top_n(n = 1) # nolint
maxf2 <- metrics_per_threshold_tsv %>% select(Threshold, `F2-score`) %>% top_n(n = 1) # nolint
melted <- melt(metrics_per_threshold_tsv, id.vars = "Threshold", 
    measure.vars = c("Precision", "Recall", "F1-score", "F2-score"));

p <- ggplot(melted, aes(x = Threshold, y = value,  colour = variable)) +
    xlim(0, 1) + ylim(0, 1);
p <- p + geom_line(size = size_geom_line) +
    scale_colour_manual(values = colours_nice) +
    standard_style +
    labs(x = opt$xname, y = opt$yname);
p <- p + geom_vline(xintercept = c(maxf1$Threshold,maxf2$Threshold), size = 1,
    colour = c("#7570b3", "#e6ab02"), linetype = "longdash")
p <- p + geom_text(
    aes(
        x = maxf1$Threshold, label = round(maxf1$Threshold, digits = 3),
        y = 0.97
    ), colour = "#7570b3", angle = 90, vjust = +1.1, size = 10)
p <- p + geom_text(
    aes(
        x = maxf2$Threshold, label = round(maxf2$Threshold, digits = 3),
        y = 0.97
    ), colour = "#e6ab02", angle = 90, vjust = -0.1, size = 10)

ggsave(p, file = opt$output, width = 10.0, height = 10)
