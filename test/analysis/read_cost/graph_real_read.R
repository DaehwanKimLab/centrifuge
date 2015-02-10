# In order to run the script, type:
#     "Rscript graph_real_read.R analysis_real_20M.data"

# install.packages("ggplot2")
library(grid)
library(ggplot2)
library(scales)
# install.packages("xlsx")
library(xlsx)

options(scipen=12)

# my_colors = c('green', 'brown', 'yellow', 'darkorange', 'red', 'magenta', 'gray')
my_colors = c('#008000A0', '#A52A2AA0', '#FFFF00A0', '#FF8C00A0', '#0000FFA0', '#FF0000A0', '#FF00FFA0', '#808080A0')
my_colors2 = c("red", "yellow", "green", "gray", "magenta", "blue")

canonical_aligner_name <- function(aligner) {
  if (aligner == "star") {
    "STAR"
  } else if (aligner == "starx2") {
    "STARx2"
  } else if(aligner == "olego") {
    "OLego"
  } else if(aligner == "gsnap") {
    "GSNAP"
  } else if(aligner == "hisatx1") {
    "HISATx1"
  } else if(aligner == "hisatx2") {
    "HISATx2"
  } else if(aligner == "hisat") {
    "HISAT"
  } else if(aligner == "tophat2") {
    "TopHat2"
  } else {
    "NA"    
  }
}

aligners <- c("hisatx1", "olego", "gsnap", "star", "starx2", "hisat", "hisatx2", "tophat2")
aligner_names <- c("HISATx1", "OLego", "GSNAP", "STAR", "STARx2", "HISAT", "HISATx2", "TopHat2")
types <- c("0", "1", "2", "3")

args <- commandArgs(trailingOnly = TRUE)

data_fname <- args[1]
output_fname <- paste(data_fname, ".pdf", sep = "")
output_sep_fname <- paste(data_fname, "_sep.pdf", sep = "")
xlsx_fname <- paste(data_fname, ".xls", sep = "")

previous_theme <- theme_set(theme_bw())

read_alignment <- read.table(data_fname, header = TRUE, sep = "\t")
read_alignment <- read_alignment[read_alignment$aligner != "bowtie" & read_alignment$aligner != "bowtie2" & read_alignment$edit_distance <= 3,
                                c("end_type", "aligner", "edit_distance", "mapped_reads", "junction_reads", "gtf_junction_reads", "runtime")]

# print(read_alignment)

num_rows = nrow(read_alignment)

read_overall <- data.frame(matrix(nrow = num_rows, ncol = 3,
                                    dimnames = list(NULL, c("aligner", "type", "alignment"))))
read_junction <- data.frame(matrix(nrow = num_rows, ncol = 4,
                                   dimnames = list(NULL, c("aligner", "type", "junctions", "gtf_junctions"))))

j = 1
for (i in 1:num_rows) {
  aligner <- canonical_aligner_name(as.character(read_alignment$aligner[i]))

  read_overall$aligner[i] <- aligner
  read_overall$type[i] <- read_alignment$edit_distance[i]
  read_overall$alignment[i] <- read_alignment$mapped_reads[i]
  
  read_junction$aligner[i] <- aligner
  read_junction$type[i] <- read_alignment$edit_distance[i]
  read_junction$junctions[i] <- read_alignment$junction_reads[i]
  read_junction$gtf_junctions[i] <- read_alignment$gtf_junction_reads[i]
}

read_overall$alignerf <- factor(read_overall$aligner, level = aligner_names)
read_overall$typef <- factor(read_overall$type, level = types)
# print(read_overall)

read_junction$alignerf <- factor(read_junction$aligner, level = aligner_names)
read_junction$typef <- factor(read_junction$type, level = types)
# print(read_junction)

min_num_alignments = min(read_overall[, "alignment"])
max_num_alignments = max(read_overall[, "alignment"])

limit_min_alignments = round(min_num_alignments / 5000000) * 5000000 - 5000000
limit_max_alignments = round(max_num_alignments / 5000000) * 5000000 + 5000000

min_num_spliced_alignments = min(read_junction[, "gtf_junctions"])
max_num_spliced_alignments = max(read_junction[, "gtf_junctions"])

limit_min_spliced_alignments = round(min_num_spliced_alignments / 5000000) * 5000000 - 5000000
limit_max_spliced_alignments = round(max_num_spliced_alignments / 5000000) * 5000000 + 5000000

my_plot1 <- ggplot(data = read_overall, aes(x = alignerf, y = alignment, fill = alignerf)) +
  geom_bar(stat = "identity", width = 0.4) +
  xlab("") +
  ylab("Cumulative number of alignments") +
  # scale_y_continuous(limits = c(limit_min_alignments, limit_max_alignments)) +
  scale_y_continuous(labels = comma) + 
  coord_cartesian(ylim = c(limit_min_alignments, limit_max_alignments)) +
  geom_text(aes(label = comma(alignment)), size = 4, angle = -90) + 
  facet_grid(facets = . ~ typef) +
  # scale_fill_manual(values=my_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray"), panel.border = element_blank()) +
  theme(axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none")

my_plot2 <- ggplot(data = read_junction, aes(x = alignerf, y = gtf_junctions, fill = alignerf)) +
  geom_bar(stat = "identity", width = .4) +
  xlab("") +
  ylab("Cumulative number of spliced alignments") +
  scale_y_continuous(labels = comma) + 
  coord_cartesian(ylim = c(limit_min_spliced_alignments, limit_max_spliced_alignments)) +
  geom_text(aes(label = comma(gtf_junctions)), size = 4, angle = -90) +
  facet_grid(facets = . ~ typef) +
  # scale_fill_manual(values=my_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray"), panel.border = element_blank()) +
  theme(axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none")


pdf(file = output_fname, width = 15, height = 16)
pushViewport(viewport(layout = grid.layout(2, 1, heights = unit(c(.5, .5), c("null", "null")))))
print(my_plot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(my_plot2, vp = viewport(layout.pos.row = 2, layout.pos.col = 1))
popViewport()
dev.off()

pdf(file = output_sep_fname, width = 15, height = 10)
print(my_plot1)
print(my_plot2)
dev.off()

write.xlsx(read_overall[, c("aligner", "type", "alignment")], xlsx_fname, sheetName = "read")
write.xlsx(read_junction[, c("aligner", "type", "gtf_junctions")], xlsx_fname, sheetName = "junction read", append = TRUE)


theme_set(previous_theme)
