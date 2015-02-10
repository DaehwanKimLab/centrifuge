# In order to run the script, type:
#     "Rscript graph_read_cost.R genome_sim_20M_single.analysis"

# install.packages("ggplot2")
library(grid)
library(ggplot2)
library(scales)
# install.packages("xlsx")
library(xlsx)

# my_colors = c('green', 'brown', 'yellow', 'darkorange', 'red', 'magenta', 'gray')
my_colors = c('#008000A0', '#A52A2AA0', '#FFFF00A0', '#FF8C00A0', '#0000FFA0', '#FF0000A0', '#FF00FFA0', '#808080A0')
my_colors2 = c("red", "yellow", "green", "gray", "magenta", "blue")
pie_colors = c("red", "turquoise", "orange", "blue", "purple")

canonical_aligner_name <- function(aligner) {
  if (aligner == "star") {
    "STAR"
  } else if (aligner == "starx2") {
    "STARx2"
  } else if(aligner == "olego") {
    "OLego"
  } else if(aligner == "gsnap") {
    "GSNAP"
  } else if(aligner == "hisat") {
    "HISAT"
  } else if(aligner == "hisatx2") {
    "HISATx2"
  } else if(aligner == "hisatx1") {
    "HISATx1"
  } else if(aligner == "tophat2") {
    "TopHat2"
  } else if(aligner == "tophat2_0") {
    "TopHat2"
  } else {
    "NA"    
  }
}


aligners <- c("hisatx1", "olego", "gsnap", "star", "starx2", "hisat", "hisatx2", "tophat2")
aligner_names <- c("HISATx1", "OLego", "GSNAP", "STAR", "STARx2", "HISAT", "HISATx2", "TopHat2")
types <- c("all", "M", "2M_gt_15", "2M_8_15", "2M_1_7", "gt_2M")
mapped_types <- c("Correctly and uniquely mapped", "Correctly mapped (multimapped)", "Incorrectly mapped", "Unmapped")

args <- commandArgs(trailingOnly = TRUE)

data_fname <- args[1]
output_fname <- paste(data_fname, ".pdf", sep = "")
xlsx_fname <- paste(data_fname, ".xls", sep = "")

read_cost <- read.table(data_fname, header = TRUE, sep = "\t")
num_total_reads <- unique(read_cost[read_cost$type == "all", "num_reads"])

pdf(file = output_fname, width = 15, height = 10)

previous_theme <- theme_set(theme_bw())
reads <- read_cost[read_cost$type != "all", c("type", "num_reads")]
reads <- unique(reads)

# num_etc_reads <- num_total_reads - sum(reads[,2])
# etc_reads <- data.frame(type = "etc", all = num_etc_reads)
# reads <- rbind(reads, etc_reads)
  
#layout(matrix(c(1,1,2,1,1,3,0,0,4), 3, 3, byrow = TRUE))
#layout(matrix(c(0,1,1,0,0,1,1,0,2,3,4,5), 3, 4, byrow = TRUE))
layout(matrix(c(1,1,1,2,3, 1,1,1,4,5), 2, 5, byrow = TRUE))

# http://www.statmethods.net/graphs/pie.html
slices <- reads[,2] / num_total_reads
lbls <- round(slices * 1000) / 10
lbls <- paste(lbls, "%", sep="") # ad % to labels
lbls <- paste(lbls, "(")
lbls <- paste(lbls, reads[,1], sep="")
lbls <- paste(lbls, ")", sep="")
# pie(slices, labels = lbls, main = "Reads", col = rainbow(length(lbls)), cex=1.3)
pie(slices, labels = lbls, main = "Reads", col = pie_colors, cex=1.3)

write.xlsx(reads, xlsx_fname, sheetName = "reads") 

pie_runtime <- read_cost[read_cost$aligner != "bowtie" &
                         read_cost$aligner != "bowtie2" &
                         read_cost$type != "all",
                         c("type", "aligner", "time")]

# print(pie_runtime)

for (aligner_idx in 1:length(aligners)) {
  aligner <- aligners[aligner_idx]
  if(aligner == "hisat" | aligner == "hisatx2" | aligner == "tophat2" | aligner == "starx2") {
    next
  }
  
  aligner_runtime <- pie_runtime[pie_runtime$aligner == aligner, c("type", "time")]
    
  # print(aligner_runtime)
    
  slices <- aligner_runtime[,2] / sum(aligner_runtime[,2])
  lbls <- round(slices * 1000) / 10
  lbls <- paste(lbls, "%", sep="") # ad % to labels
  # lbls <- paste(lbls, "(")
  # lbls <- paste(lbls, aligner_runtime[,1], sep="")
  #lbls <- paste(lbls, ")", sep="")

  aligner_name <- canonical_aligner_name(aligner)
        
  # pie(slices, labels = lbls, main = paste("Relative runtime by", aligner_name), col = rainbow(length(lbls)), cex=1.3)
  pie(slices, labels = lbls, main = paste("Relative runtime by", aligner_name), col = pie_colors, cex=1.3)
}

temp_read_cost <- read_cost
read_runtime <- temp_read_cost[temp_read_cost$aligner != "bowtie" & temp_read_cost$aligner != "bowtie2",
                               c("end_type", "type", "aligner", "num_reads", "time", "mapped_reads", "unique_mapped_reads", "incorrect_mapped_reads", "unmapped_reads")]

read_runtime$multi_mapped_reads <- read_runtime$mapped_reads - read_runtime$unique_mapped_reads

end_type = "Reads"
if(read_runtime$end_type[1] == "paired") {
   end_type = "Pairs"
}

num_rows = nrow(read_runtime)

read_alignment <- data.frame(matrix(nrow = num_rows, ncol = 3,
                                    dimnames = list(NULL, c("aligner", "type", "alignment"))))
read_alignment2 <- data.frame(matrix(nrow = num_rows * 4, ncol = 7,
               	                     dimnames = list(NULL, c("aligner", "type", "mapped_type", "alignment", "alignment_multi", "num_read", "num_alignment"))))
read_runtime2 <- data.frame(matrix(nrow = num_rows - 5 * 4, ncol = 6,
                                   dimnames = list(NULL, c("aligner", "type", "type2", "num_read", "time", "read_per_sec"))))

j = 1
for (i in 1:num_rows) {
  aligner <- canonical_aligner_name(as.character(read_runtime$aligner[i]))
  type <- as.character(read_runtime$type[i])

  read_alignment$aligner[i] <- aligner
  read_alignment$type[i] <- type
  read_alignment$alignment[i] <- round(read_runtime$mapped_reads[i] / read_runtime$num_reads[i] * 1000) / 10

  for (i2 in 1:4) {
    read_alignment2$aligner[(i-1)*4+i2] <- aligner
    read_alignment2$type[(i-1)*4+i2] <- type
  }

  read_alignment2$mapped_type[(i-1)*4+1] <- "Correctly and uniquely mapped"
  read_alignment2$alignment[(i-1)*4+1] <- round(read_runtime$unique_mapped_reads[i] / read_runtime$num_reads[i] * 1000) / 10
  read_alignment2$alignment_multi[(i-1)*4+1] <- round(read_runtime$mapped_reads[i] / read_runtime$num_reads[i] * 1000) / 10
  read_alignment2$num_read[(i-1)*4+1] <- read_runtime$num_reads[i]
  read_alignment2$num_alignment[(i-1)*4+1] <- read_runtime$unique_mapped_reads[i]
  read_alignment2$mapped_type[(i-1)*4+2] <- "Correctly mapped (multimapped)"
  read_alignment2$alignment[(i-1)*4+2] <- round((read_runtime$mapped_reads[i] - read_runtime$unique_mapped_reads[i]) / read_runtime$num_reads[i] * 1000) / 10
  read_alignment2$num_read[(i-1)*4+2] <- read_runtime$num_reads[i]
  read_alignment2$num_alignment[(i-1)*4+2] <- read_runtime$mapped_reads[i] - read_runtime$unique_mapped_reads[i]
  read_alignment2$mapped_type[(i-1)*4+3] <- "Incorrectly mapped"
  read_alignment2$alignment[(i-1)*4+3] <- round(read_runtime$incorrect_mapped_reads[i] / read_runtime$num_reads[i] * 1000) / 10
  read_alignment2$num_read[(i-1)*4+3] <- read_runtime$num_reads[i]
  read_alignment2$num_alignment[(i-1)*4+3] <- read_runtime$incorrect_mapped_reads[i]
  read_alignment2$mapped_type[(i-1)*4+4] <- "Unmapped"
  read_alignment2$alignment[(i-1)*4+4] <- round(read_runtime$unmapped_reads[i] / read_runtime$num_reads[i] * 1000) / 10
  read_alignment2$num_read[(i-1)*4+4] <- read_runtime$num_reads[i]
  read_alignment2$num_alignment[(i-1)*4+4] <- read_runtime$unmapped_reads[i]

  if(type == "all" | (aligner != "STARx2" & aligner != "HISAT" & aligner != "HISATx2" & aligner != "TopHat2")) {
    read_runtime2$aligner[j] <- aligner
    read_runtime2$type[j] <- type
    read_runtime2$type2[j] <- "all"
    read_runtime2$num_read[j] <- read_runtime$num_reads[i]
    read_runtime2$time[j] <- read_runtime$time[i]
    read_runtime2$read_per_sec[j] <-
      min(read_runtime2$num_read[j], round(read_runtime2$num_read[j] / read_runtime2$time[j]))
    j <- j + 1
  }
}

for (i in 1:length(types)) {
  read_runtime2[read_runtime2$type == types[i], "read_per_sec_annotated"] <- 
    round(max(read_runtime2[read_runtime2$type == types[i], "read_per_sec"]) / 2)
}
    
read_alignment$alignerf <- factor(read_alignment$aligner, level = aligner_names)
read_alignment$typef <- factor(read_alignment$type, level = types)
# print(read_alignment)

read_alignment2$alignerf <- factor(read_alignment2$aligner, level = aligner_names)
read_alignment2$typef <- factor(read_alignment2$type, level = types)
read_alignment2$mapped_typef <- factor(read_alignment2$mapped_type, level = mapped_types)
# print(read_alignment2)

read_runtime2$alignerf <- factor(read_runtime2$aligner, level = aligner_names)
read_runtime2$typef <- factor(read_runtime2$type, level = types)
# print(read_runtime2)

for (aligner_idx in 1:length(aligners)) {
  aligner <- aligners[aligner_idx]
}

ylab_text = paste(end_type, "processed per second")

my_plot1 <- ggplot(data = read_runtime2[read_runtime2$type2 == "all",], aes(x = alignerf, y = read_per_sec, fill = alignerf)) +
  geom_bar(stat = "identity", width = .3) +
  xlab("") +
  ylab(ylab_text) + 
  geom_text(aes(y = read_per_sec_annotated, label = comma(read_per_sec)), size = 4) +
  # facet_grid(facets = typef ~ ., scales = "free_y") +
  facet_grid(facets = typef ~ .) +
  # scale_fill_manual(values=my_colors) +
  scale_y_continuous(labels = comma) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray"), panel.border = element_blank()) + 
  theme(text=element_text(size=15), axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none")

my_plot2 <- ggplot(data = read_alignment2[,], aes(x = alignerf, y = alignment, fill = mapped_typef)) +
  geom_bar(stat = "identity", width = 0.3) +
  xlab("") +
  # ylab("Alignment sensitivity (%)") +
  ylab("") + 
  geom_text(aes(y = 50.0, label = ifelse(mapped_typef == "Correctly and uniquely mapped", alignment, "")), size = 4) +
  geom_text(aes(y = 36.0, label = ifelse(mapped_typef == "Correctly and uniquely mapped", paste("(", alignment_multi, ")", sep = ""), "")), size = 3.2) +
  facet_grid(facets = typef ~ .) +
  scale_fill_manual(values=my_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray"), panel.border = element_blank()) + 
  theme(text=element_text(size=15), axis.text.x=element_text(angle=-90)) +
  theme(legend.position = c(0.5, -0.17), legend.title = element_blank())
  # theme(label.text=element_text(size=10))


grid.newpage()
pushViewport(
  viewport(layout = grid.layout(
  		      2, 2, widths = unit(c(.52, .48), c("null", "null")), heights = unit(c(.9, .1), c("null", "null"))
		    )
  )
)


print(my_plot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(my_plot2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

write.xlsx(read_runtime2[read_runtime2$type2 == "all",c("aligner", "type", "num_read", "time")], xlsx_fname, sheetName = "runtime", append = TRUE)
write.xlsx(read_runtime[, c("aligner", "type", "num_reads", "unique_mapped_reads", "multi_mapped_reads", "incorrect_mapped_reads", "unmapped_reads")], xlsx_fname, sheetName = "alignment", append = TRUE)

# summary(my_plot2)
popViewport()

grid.newpage()
pushViewport(
  viewport(layout = grid.layout(
  		      2, 2, widths = unit(c(.52, .48), c("null", "null")), heights = unit(c(.4, .6), c("null", "null"))
		    )
  )
)

my_plot1_all <- ggplot(data = read_runtime2[read_runtime2$type == "all" & read_runtime2$type2 == "all",], aes(x = alignerf, y = read_per_sec, fill = alignerf)) +
  geom_bar(stat = "identity", width = .3) +
  xlab("") +
  ylab(ylab_text) +
  # daehwan - for debugging purposes
  # geom_text(aes(y = read_per_sec_annotated, label = comma(read_per_sec)), size = 5) +
  geom_text(aes(y = read_per_sec + 8000, label = comma(read_per_sec)), size = 5) +
  # facet_grid(facets = typef ~ ., scales = "free_y") +
  # facet_grid(facets = typef ~ .) +
  # scale_fill_manual(values=my_colors) +
  scale_y_continuous(labels = comma) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray"), panel.border = element_blank()) + 
  theme(text=element_text(size=15), axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none")

print(my_plot1_all, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

popViewport()

write.xlsx(read_runtime2[read_runtime2$type == "all",c("aligner", "num_read", "time")], xlsx_fname, sheetName = "runtime_all", append = TRUE)

grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(.68, .32), c("null", "null")))))

print(my_plot2, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

popViewport()


# plot mapped_reads and unique_mapped_reads
grid.newpage()
pushViewport(
  viewport(layout = grid.layout(
		      2, 2, widths = unit(c(.54, .46), c("null", "null")), heights = unit(c(.5, .5), c("null", "null"))
		    )
  )
)

my_plot2_all <- ggplot(data = read_alignment2[read_alignment2$type == "all",], aes(x = alignerf, y = alignment, fill = mapped_typef)) +
  geom_bar(stat = "identity", width = 0.3) +
  xlab("") +
  # ylab("Alignment sensitivity (%)") +
  ylab("") +
  # daehwan - for debugging purposes
  # geom_text(aes(y = 50.0, label = ifelse(mapped_typef == "Correctly and uniquely mapped", alignment, "")), size = 5) +
  # geom_text(aes(y = 42.0, label = ifelse(mapped_typef == "Correctly and uniquely mapped", paste("(", alignment_multi, ")", sep = ""), "")), size = 4) +
  geom_text(aes(y = 114, label = ifelse(mapped_typef == "Correctly and uniquely mapped", alignment, "")), size = 5) +
  geom_text(aes(y = 106, label = ifelse(mapped_typef == "Correctly and uniquely mapped", paste("(", alignment_multi, ")", sep = ""), "")), size = 4) +
  scale_y_continuous(breaks=seq(0,100,25)) + 
  # facet_grid(facets = typef ~ .) +
  scale_fill_manual(values=my_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray"), panel.border = element_blank()) + 
  theme(text=element_text(size=15), axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "top", legend.title=element_blank())

print(my_plot2_all, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

popViewport()

write.xlsx(read_runtime[read_runtime$type == "all", c("aligner", "num_reads", "unique_mapped_reads", "multi_mapped_reads", "incorrect_mapped_reads", "unmapped_reads")], xlsx_fname, sheetName = "alignment_all", append = TRUE)


grid.newpage()
pushViewport(
  viewport(layout = grid.layout(
		      2, 2, widths = unit(c(.58, .42), c("null", "null")), heights = unit(c(.7, .3), c("null", "null"))
		    )
  )
)

my_plot2_all <- ggplot(data = read_alignment2[read_alignment2$type == "2M_8_15" | read_alignment2$type == "2M_1_7",], aes(x = alignerf, y = alignment, fill = mapped_typef)) +
  geom_bar(stat = "identity", width = 0.3) +
  xlab("") +
  # ylab("Alignment sensitivity (%)") +
  ylab("") +
  # daehwan - for debugging purposes
  # geom_text(aes(y = 50.0, label = ifelse(mapped_typef == "Correctly and uniquely mapped", alignment, "")), size = 5) +
  # geom_text(aes(y = 39.0, label = ifelse(mapped_typef == "Correctly and uniquely mapped", paste("(", alignment_multi, ")", sep = ""), "")), size = 4) +
  geom_text(aes(y = 114.0, label = ifelse(mapped_typef == "Correctly and uniquely mapped", alignment, "")), size = 5) +
  geom_text(aes(y = 106.0, label = ifelse(mapped_typef == "Correctly and uniquely mapped", paste("(", alignment_multi, ")", sep = ""), "")), size = 4) +
  scale_y_continuous(breaks=seq(0,100,25)) + 
  facet_grid(facets = typef ~ .) +
  scale_fill_manual(values=my_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), axis.line = element_line(colour = "gray"), panel.border = element_blank()) + 
  theme(text=element_text(size=15), axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "top", legend.title = element_blank())

print(my_plot2_all, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

popViewport()

write.xlsx(read_runtime[read_runtime$type == "2M_8_15" | read_runtime$type == "2M_1_7", c("aligner", "type", "num_reads", "unique_mapped_reads", "multi_mapped_reads", "incorrect_mapped_reads", "unmapped_reads")], xlsx_fname, sheetName = "alignment_small", append = TRUE)


# draw a graph for sensitivity and precision of splice sites
grid.newpage()
pushViewport(viewport(layout = grid.layout(1, 2, widths = unit(c(.5, .5), c("null", "null")))))

temp_read_cost <- read_cost
read_junction <- temp_read_cost[temp_read_cost$aligner != "bowtie" & temp_read_cost$aligner != "bowtie2" & temp_read_cost$type != "M",
                                c("type", "aligner", "true_gtf_junctions", "temp_junctions", "temp_gtf_junctions")]

num_rows = nrow(read_junction)

read_junction2 <- data.frame(matrix(nrow = num_rows - 3 * 4, ncol = 4,
                             dimnames = list(NULL, c("aligner", "type", "sensitivity", "precision"))))

j = 1
for (i in 1:num_rows) {
  aligner <- canonical_aligner_name(as.character(read_junction$aligner[i]))
  type <- as.character(read_junction$type[i])

  if(type == "all" | (aligner != "HISAT" & aligner != "HISATx2" & aligner != "TopHat2")) {
    read_junction2$aligner[j] <- aligner
    read_junction2$type[j] <- type
    read_junction2$sensitivity[j] <- round(read_junction$temp_gtf_junctions[i] / read_junction$true_gtf_junctions[i] * 1000) / 10
    read_junction2$accuracy[j] <- round(read_junction$temp_gtf_junctions[i] / read_junction$temp_junctions[i] * 1000) / 10
    j <- j + 1
  }
}

read_junction2$alignerf <- factor(read_junction2$aligner, level = aligner_names)
read_junction2$typef <- factor(read_junction2$type, level = types)

my_plot1 <- ggplot(data = read_junction2, aes(x = alignerf, y = sensitivity, fill = alignerf)) +
  geom_bar(stat = "identity", width = 0.3) +
  xlab("") +
  ylab("Sensitivity of Splice Sites (%)") + 
  geom_text(aes(y = 50.0, label = sensitivity), size = 5) + 
  facet_grid(facets = typef ~ .) +
  scale_fill_manual(values=my_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text=element_text(size=15), axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none")

print(my_plot1, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))

my_plot2 <- ggplot(data = read_junction2, aes(x = alignerf, y = accuracy, fill = alignerf)) +
  geom_bar(stat = "identity", width = 0.3) +
  xlab("") +
  ylab("Accuracy of Splice Sites (%)") + 
  geom_text(aes(y = 50.0, label = accuracy), size = 5) + 
  facet_grid(facets = typef ~ .) +
  scale_fill_manual(values=my_colors) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank()) + 
  theme(text=element_text(size=15), axis.text.x=element_text(angle=-90)) +
  theme(legend.position = "none")

print(my_plot2, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))

dev.off()

theme_set(previous_theme)
