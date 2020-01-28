#!/usr/bin/env Rscript

# Thanks to Matthew Hall, the code's original author.

list.of.packages <- c("ggplot2", "reshape", "grid", "gridExtra", "gtable", "scales", "argparse", "assertthat")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)){
  cat("Missing dependencies; replacing the path below appropriately, run\nsudo [path to your shiver code]/tools/AlignmentPlotting_PackageInstall.R\nthen try again.\n")
  quit(save="no", status=1)
}


suppressMessages(require(ggplot2, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(reshape, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(grid, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gridExtra, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(gtable, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(scales, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(argparse, quietly=TRUE, warn.conflicts=FALSE))
suppressMessages(require(assertthat, quietly=TRUE, warn.conflicts=FALSE))



# Miseq
gene.height <- 3.2
other.seq.height = 3.5
gene.font.size <- 2.7
v.loop.font.size <- 2.5
plot.relative.heights <- c(2,4,6.5)
plot.width <- 13
plot.total.height <- 3.8
legend.position <- "none"

# Hiseq
#gene.height <- 3.2
#other.seq.height <- 3.6
#gene.font.size <- 2.7
#v.loop.font.size <- 2.5
#plot.relative.heights <- c(2,5,6.5)
#plot.width <- 13
#plot.total.height <- 5
#legend.position <- "none"


arg_parser = ArgumentParser(description=paste("Use this script to produce a",
"plot showing an alignment of sequences, genes relative to those sequences,",
"and the coverage obtained by mapping to two of those sequences (the alignment",
"may contain additional sequences with those two). NOTE: I am not",
"able to offer support for this script as I am for all of the other code in",
"shiver. It was written, not by me, with the intention of one-off use for the",
"figures in the shiver paper. In my limited experience of R code I find it",
"produces incomprehensible error messages; I wish you luck debugging this code",
"if it doesn't work for you. Also note that there are lots of plotting",
"parameters - sizes, placements, colours etc. - you can modify these by",
"opening this file and modifying the code."))

arg_parser$add_argument("coverageFile",
help=paste("A csv file of the format output by shiver's",
"tools/AlignBaseFreqFiles.py script with the option --coverage-only. i.e. the",
"first column should be position in the alignment,",
"and the fourth and fifth columns should be coverage with respect",
"to the two references. (The first line of the csv should be the field",
"names)."))  
arg_parser$add_argument("coloursFile", help="A csv file produced by running shiver's tools/ConvertAlnToColourCodes.py on your alignment.")  
arg_parser$add_argument("GeneCoordsFile", help='A csv file with the following fields: "ReadingFrame,GeneName,StartPos,EndPos". The start and end positions should be in the coordinates of your alignment. An optional extra field "OnTop" should contain binary values "yes" or "no", specifying whether that gene should be displayed on top of another gene at the same position (we use yes values for putting loop regions of the env gene on top of env itself). If your alignment contains HXB2, you could create this file by running a command like this (assuming your shiver code lives in ~/shiver/): echo "ReadingFrame,GeneName,StartPos,EndPos" > MyGeneCoordsFile.csv; while read line; do read name RF seq <<< $(echo $line); start=$(~/shiver/tools/FindSubSeqsInAlignment.py "$i" B.FR.83.HXB2_LAI_IIIB_BRU.K03455 --start "$seq" --alignment-coords); end=$(~/shiver/tools/FindSubSeqsInAlignment.py "$i" B.FR.83.HXB2_LAI_IIIB_BRU.K03455 --end "$seq" --alignment-coords); echo "$RF,$name,$start,$end"; done < ~/shiver/info/HXB2_GeneName_ReadingFrame_Sequence.txt >> MyGeneCoordsFile.csv')  
arg_parser$add_argument("outputFile", help='(Will be in pdf format.)')  

args <- arg_parser$parse_args()

cov.file.name <- args$coverageFile
col.file.name <- args$coloursFile
gen.file.name <- args$GeneCoordsFile
out.file.name <- args$outputFile

AlignPlots <- function(...) {
  LegendWidth <- function(x) x$grobs[[8]]$grobs[[1]]$widths[[4]]
  
  plots.grobs <- lapply(list(...), ggplotGrob)
  
  max.widths <- do.call(unit.pmax, lapply(plots.grobs, "[[", "widths"))
  plots.grobs.eq.widths <- lapply(plots.grobs, function(x) {
    x$widths <- max.widths
    x
  })
  
  legends.widths <- lapply(plots.grobs, LegendWidth)
  max.legends.width <- do.call(max, legends.widths)
  plots.grobs.eq.widths.aligned <- lapply(plots.grobs.eq.widths, function(x) {
    if (is.gtable(x$grobs[[8]])) {
      x$grobs[[8]] <- gtable_add_cols(x$grobs[[8]],
                                      unit(abs(diff(c(LegendWidth(x),
                                                      max.legends.width))),
                                           "mm"))
    }
    x
  })
  
  plots.grobs.eq.widths.aligned
}

process.string <- function(string){
  anything.exists <- F
  
  gs <- gregexpr("g*", string)
  lengths <- attr(gs[[1]], "match.length")
  nonzero.lengths <- which(lengths!=0)
  if(length(nonzero.lengths)>0){
    starts <- gs[[1]][nonzero.lengths]
    display.starts <- starts - 0.5
    ends <- starts + lengths[nonzero.lengths] - 1
    display.ends <- ends + 0.5
    out.g <- data.frame(start = starts, end = ends, display.start = display.starts, display.end = display.ends, char = "g")
    
    anything.exists <- T
    out <- out.g
  }
  
  ds <- gregexpr("d*", string)
  lengths <- attr(ds[[1]], "match.length")
  nonzero.lengths <- which(lengths!=0)
  if(length(nonzero.lengths)>0){
    starts <- ds[[1]][nonzero.lengths]
    display.starts <- starts - 0.5
    ends <- starts + lengths[nonzero.lengths] - 1
    display.ends <- ends + 0.5
    out.d <- data.frame(start = starts, end = ends, display.start = display.starts, display.end = display.ends, char = "d")
    
    if(anything.exists){
      out <- rbind(out, out.d)
    } else {
      anything.exists <- T
      out <- out.d
    }
  }
  
  bs <- gregexpr("b*", string)
  lengths <- attr(bs[[1]], "match.length")
  nonzero.lengths <- which(lengths!=0)
  if(length(nonzero.lengths)>0){
    starts <- bs[[1]][nonzero.lengths]
    display.starts <- starts - 0.5
    ends <- starts + lengths[nonzero.lengths] - 1
    display.ends <- ends + 0.5
    out.b <- data.frame(start = starts, end = ends, display.start = display.starts, display.end = display.ends, char = "b")
    
    if(anything.exists){
      out <- rbind(out, out.b)
    } else {
      anything.exists <- T
      out <- out.b
    }
  }
  
  return(out)
}

read.csv.and.check.col.names <- function(csv.file, all.required.col.names,
                                         discard.unwanted.cols=FALSE,
                                         check.names=TRUE,
                                         stringsAsFactors=FALSE) {
  
  # Read the file, check all required columns present
  assert_that(file.exists(csv.file))
  data.frame <- read.table(csv.file, sep=",", header=T, check.names=check.names,
  stringsAsFactors=stringsAsFactors)
  col.names <- colnames(data.frame)
  for (expected.col.name in all.required.col.names) {
    if (! expected.col.name %in% col.names) {
      stop("Expected column name ", expected.col.name, " missing from ",
           csv.file)
    }
  }
  
  if (discard.unwanted.cols) {
    data.frame <- data.frame[all.required.col.names]
  }
  
  return(data.frame)
}


cov.data <- read.table(cov.file.name, sep=",", header=T, stringsAsFactors = F)
cov.data <- cov.data[,c(1,4,5)]
comp.factors <- c(colnames(cov.data)[c(3,2)])
cov.data <- melt(cov.data, id = 1)
cov.data$variable <- factor(cov.data$variable, comp.factors)

pos.data <- read.table(col.file.name, sep=",", header=F, stringsAsFactors = F)
colnames(pos.data) <- c("Name", "Sequence")
pos.data$Name <- factor(pos.data$Name, rev(pos.data$Name))

line.list <- lapply(pos.data$Sequence, process.string)

first <- T

for(sequence.no in seq(1, nrow(pos.data))){
  temp <- line.list[[sequence.no]]
  temp$phantom.pos <- sequence.no
  temp$sequence.name <- pos.data$Name[sequence.no]
  if(first){
    line.df <- temp
    first <- F
  } else {
    line.df <- rbind(line.df, temp)
  }
  
}

ann.data <- read.csv.and.check.col.names(gen.file.name, c('ReadingFrame','GeneName','StartPos','EndPos'))
ann.data$display.start <- ann.data$StartPos-0.5
ann.data$display.end <- ann.data$EndPos+0.5
ann.data <- ann.data[order(ann.data$GeneName),]


stack.genes <- "OnTop" %in% colnames(ann.data)
if (stack.genes) {
  ann.data.main <- ann.data[which(ann.data[,5]=="no"),]
  ann.data.extra <- ann.data[which(ann.data[,5]=="yes"),]
} else {
  ann.data.main <- ann.data
}

graph.cov <- ggplot(data=cov.data)

graph.cov <- graph.cov + 
  geom_line(aes(x=Alignment.position, y=value, colour=variable, alpha=0.5)) +
  ylab("Coverage") +
  xlab("Alignment position") +
  scale_x_continuous(limits=c(min(cov.data$Alignment.position)-1, max(cov.data$Alignment.position)+1), expand=c(0,100)) +
  scale_y_log10(limits = c(1,NA), labels = comma) +
  scale_colour_manual(values=c("red", "blue")) +
  theme_bw() +
  #theme(legend.title = element_blank()) +
  theme(plot.margin = unit(c(0,0.5,0.5,0.5), "cm"),
  legend.position=legend.position)
  #legend.text = element_text(size=10))

graph.ann <- ggplot()

graph.ann <- graph.ann + 
  geom_segment(data = ann.data.main, aes(y=ReadingFrame, yend = ReadingFrame, x = display.start, xend=display.end, colour=GeneName), size=gene.height) +
  geom_text(data = ann.data.main, aes(label=GeneName, y=ReadingFrame, x=display.start), hjust="left", nudge_x=10, size=gene.font.size) +
  theme_bw() +
  scale_x_continuous(limits=c(min(cov.data$Alignment.position)-1, max(cov.data$Alignment.position)+1), expand=c(0,100)) +
  scale_y_continuous(limits=c(0, max(ann.data$ReadingFrame)+1)) +
  scale_color_discrete(h = c(0, 720) + 15, direction = -1) +
  theme(axis.text=element_blank(),
        axis.ticks=element_blank(),
        axis.title=element_blank(),
        legend.position="none",
        panel.background=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(),
        plot.background=element_blank(),
        plot.margin=unit(c(0.5,0.5,0,0.5), "cm"),
        panel.margin=unit(c(0.5,0.5,0,0.5), "cm"))

if (stack.genes) {
graph.ann <- graph.ann +
  geom_segment(data = ann.data.extra, aes(y=ReadingFrame, yend = ReadingFrame, x = display.start, xend=display.end), size=gene.height, alpha=0.5, colour="grey10") +
  geom_text(data = ann.data.extra, aes(label=GeneName, y=ReadingFrame, x=display.start), hjust="left", nudge_x=10, size=v.loop.font.size, colour="white") 
}

graph.pos <- ggplot(data = line.df)

graph.pos <- graph.pos +
  geom_segment(aes(y=sequence.name, yend = sequence.name, x=display.start, xend = display.end, colour=char, size=char)) +
  theme_bw() +
  scale_size_manual(values=c(other.seq.height,0.5,other.seq.height)) +
  scale_color_manual(values=c("grey85", "black", "black")) +
  scale_x_continuous(limits=c(min(cov.data$Alignment.position)-1, max(cov.data$Alignment.position)+1), expand=c(0,100)) +
  theme(legend.position="none", 
        axis.ticks.y=element_blank(), 
        axis.title.y=element_blank(),
        axis.title.x=element_blank(),
        axis.ticks.x=element_blank(),
        axis.text.x=element_blank(),
        panel.grid.major.y=element_blank(),
        panel.grid.minor.y=element_blank(), 
        plot.margin=unit(c(0,0.5,0,0.5), "cm"),
        panel.margin=unit(c(0,0.5,0,0.5), "cm"))

plots1 <- AlignPlots(graph.ann, graph.pos, graph.cov)
plots1$ncol <- 1
plots1$heights <- unit(plot.relative.heights, "null")

# Append '.pdf' to the output file name if not already there.
if (! grepl('.pdf$', out.file.name)) out.file.name <- paste0(out.file.name, '.pdf')

pdf(file=out.file.name, width=plot.width, height=plot.total.height)
do.call(grid.arrange, plots1)
dev.off()
