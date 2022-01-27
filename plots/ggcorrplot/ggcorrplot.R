library(ggcorrplot)
library(tidyverse)
library(optparse)

option_list <- list(
    make_option(c("-i", "--input"), type = "character",
        help = "Input table. Multiple tables sebarated by comma.",
        dest = "input"),
    make_option(c("-c", "--column"), type = "character",
        help = "Columns for correlations, separated by comma. If not set all are used."),
    make_option(c("-b", "--bind"), action = "store_true", default = FALSE,
        help = "Bind colum."),
    make_option(c("-a", "--arrange"), type = "character",
        help = "Column for sorting. When not set no Sorting is done."),
    make_option(c("-m", "--method"), type = "character",
        help = "Correlation method. pearson or spearman (Optional default [%default])",# nolint
        default = "spearman"),
    make_option(c("-o", "--output"), type = "character",
        help = "Output file of the correlation plot")
)

arguments <- parse_args(OptionParser(option_list = option_list), positional_arguments = TRUE) # nolint

opt <- arguments$options

if (!"input" %in% names(opt)) {
  stop("--input parameter must be provided. See script usage (--help)")
}
if (!"output" %in% names(opt)) {
  stop("--output parameter must be provided. See script usage (--help)")
}


input_a <- data.frame()
input_b <- data.frame()

i <- 1
for (file in unlist(str_split(opt$input, ","))) {
    table <- read.delim(file, na.strings = "NaN")
    if ("arrange" %in% names(opt)) {
        table <- table %>% arrange(across(unlist(opt$arrange)))
    }
    if ("column" %in% names(opt)) {
        table <- table %>% select(unlist(str_split(opt$column, ",")))
    }

    if (opt$bind) {
        input_a <- rbind(input_a, table)
    } else if (i == 1) {
       input_a <- table
    } else if (i == 2) {
       input_b <- table
    } else if (i >= 2) {
       stop("More than two input files cann only be used with option --bind")
    }
    i <- i + 1
}

# these attributes are not used!
input_a <- input_a %>% select(where(is.numeric))

if (nrow(input_b) == 0) {
    correlation <- cor(input_a, method = opt$method)
} else {
    input_b <- input_b %>% select(where(is.numeric))
    correlation <- cor(input_a, input_b, method = opt$method)
}

colors <- c("#7570b3", "white", "#d95f02")

p <- ggcorrplot(correlation, method = "square",
        hc.order = TRUE,   colors = colors) +
        theme(axis.text = element_text(colour = "black"))

ggsave(p, file = opt$output, width = 10.0, height = 10.0)
