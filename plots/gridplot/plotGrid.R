library(optparse)
library(tidyverse)
library(gridExtra)



option_list <- list(
    make_option(c("-i", "--input"), type="character",
        help="Input of hyperopt result file"),
    make_option(c("-y", "--value"), type="character",
        help="Variable as value (y)."),
    make_option(c("-x", "--variable"), type="character",
        help="Variable as variable (x)."),
    make_option(c("-f", "--fill"), type="character",
        help="Variable for fill"),
    make_option(c("-g", "--group"), type="character",
        help="Variable for group"),
    make_option(c("-k", "--colour"), type="character",
        help="Variable for colour"),
    make_option(c("-c", "--cols"), type="character",
        help="Variable for cols"),
    make_option(c("-r", "--rows"), type="character",
        help="Variable for row"),
    make_option(c("-p", "--plot"), type="character",
        help="Type of plot, Only geom_boxplot, geom_path and geom_col are supported"),
    make_option(c("-o", "--output"), type="character",
        help="Output of the plot")
)
parser <- OptionParser(option_list=option_list)
arguments <- parse_args(parser, positional_arguments=TRUE)
opt <- arguments$options

if (!"input" %in% names(opt)) {
  stop("--input parameter must be provided. See script usage (--help)")
}
if (!"fill" %in% names(opt)) {
  opt$fill <- NULL
}
if (!"group" %in% names(opt)) {
  opt$group <- NULL
}
if (!"colour" %in% names(opt)) {
  opt$colour <- NULL
}
if (!"value" %in% names(opt)) {
  stop("--value parameter must be provided. See script usage (--help)")
}
if (!"output" %in% names(opt)) {
  stop("--output parameter must be provided. See script usage (--help)")
}

data <- read.delim(opt$input,stringsAsFactors=TRUE, header=TRUE) #%>% replace_na("unknown")

data[,opt$fill] <- as.factor(data[,opt$fill])

getFacetGrid <- function(){
    if ("cols" %in% names(opt) && "rows" %in% names(opt)){
        return(facet_grid(cols = vars(data[[opt$cols]]), rows= vars(data[[opt$rows]])))
    } else if ("cols" %in% names(opt)) {
        return(facet_grid(cols = vars(data[[opt$cols]])))
    } else if ("rows" %in% names(opt)) {
        return(facet_grid(cols = vars(data[[opt$rows]])))
    }
}



getPlot <- function(plot) {
    result = switch(
        plot,
        "geom_boxplot"=geom_boxplot(aes_string(y=opt$value,fill=opt$fill)),
        "geom_col"=geom_col(aes_string(x=opt$variable,y=opt$value,fill=opt$fill,group=opt$group)),
        "geom_path"=geom_path(aes_string(x=opt$variable,y=opt$value,group=opt$group,fill=opt$fill,colour=opt$colour)),
    )
    return(result)
}

p <- data %>%
    ggplot() +
    getPlot(opt$plot) +
    ylim(min(data[[opt$value]]),max(data[[opt$value]])) +
    theme_bw() + theme(legend.position="top") + 
    getFacetGrid()

ggsave(opt$output, plot=p, height=20,width=20)
