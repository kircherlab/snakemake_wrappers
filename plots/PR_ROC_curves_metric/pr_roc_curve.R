library(tidyverse)
library(reshape2)
library(optparse)

option_list <- list(
    make_option(c("-i", "--input"), type = "character",
        help = "Input table. Multiple tables sebarated by comma.",
        dest = "input"),
    make_option(c("-n", "--name"), type = "character",
        help = "Input names. Multiple seperated by comma.",
        dest = "name"),
    make_option(c("-t", "--type"), type = "character",
        default = "ROC",
        help = "Type of the plot. ROC or PR",
        dest = "type"),
    make_option(c("-x", "--xname"), type = "character",
        default = "False positive rate",
        help = "X-lab name."),
    make_option(c("-y", "--yname"), type = "character",
        default = "True positive rate",
        help = "Y-lab dname."),
    make_option(c("-c", "--labelcolumns"), type = "integer",
        default = 1,
        help = "Number of columns for label",
        dest = "label_columns"),
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
if (!opt$type %in% c("ROC", "PR")) {
  stop("--type must be PR or ROC. See script usage (--help)")
}



colours_nice <- c("#1b9e77","#d95f02","#7570b3","#e6ab02","#e7298a","#a6761d","#000000") # nolint
size_text <- 25;
size_legend <- 25;
size_line <- 1;
size_geom_line <- 2;
standard_style <- theme_bw() + theme(axis.text = element_text(colour = "black",size=size_text), axis.title = element_text(colour = "black",size=size_text), plot.title = element_text(size = size_text, face="bold"), panel.grid.major = element_blank() , panel.grid.minor = element_blank(), panel.border = element_blank(), axis.ticks = element_line(colour = "black", size=1), axis.line = element_line(colour = "black", size=size_line), legend.key =  element_blank(), legend.text = element_text(size=size_legend), legend.position="top", legend.direction="horizontal", legend.key.size = unit(2, 'lines'), legend.title=element_blank())+ theme(axis.line.x = element_line(color="black", size = size_line), axis.line.y = element_line(color="black", size = size_line)) # nolint


i <- 1
data <- data.frame()
names <- unlist(str_split(opt$name, ","))
for (file in unlist(str_split(opt$input, ","))) {
    table <- read_delim(file, delim = "\t",
        escape_double = FALSE, trim_ws = TRUE)
    table$Name <- names[i] # nolint
    data <- rbind(data, table)
    i <- i + 1
}


switch(opt$type,
    ROC = {
        p <- ggplot(data,
            aes(x = `False-positive-rate`, y = `Recall`, colour = Name)
        ) +
        geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1),
            color = "black", linetype = "dashed", size = size_geom_line) +
            labs(x = opt$xname, y = opt$yname, colour = "Name")
    },
    PR = {
        p <- ggplot(data,
            aes(x = Recall, y = Precision, colour = Name)
        )
        # geom_segment(aes(x = 0, xend = 1, y = 0.5, yend = 0.5),
            # color = "black", linetype = "dashed", size = size_geom_line)
    }, {
        stop("--type must be PR or ROC. See script usage (--help)")
    }
)

p <- p + geom_line(size = size_geom_line) +
standard_style +
guides(col = guide_legend(nrow = ceiling(length(names) / opt$label_columns),
    ncol = opt$label_columns, byrow = FALSE)) +
scale_colour_manual(values = colours_nice)

ggsave(p, file = opt$output, width = 10.0, height = 10)

# p <- p + geom_line(size = size_geom_line) +
# standard_style +
# guides(col = guide_legend(nrow = 3, byrow = FALSE)) +
# scale_colour_manual(values = colours_nice)

# ggsave(p, file = opt$output, width = 12.0, height = 12.0)
