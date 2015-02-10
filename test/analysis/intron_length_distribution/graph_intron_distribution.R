# install.packages("ggplot2")
library(grid)
library(ggplot2)
library(scales)

intron_lengths <- read.table("intron_distribution.data", header = TRUE, sep = "\t")
print(intron_lengths)

filename = paste("intron_lengths", ".pdf", sep = "")
pdf(file = filename, width = 12, height = 8)

previous_theme <- theme_set(theme_bw())

cum_count = cumsum(intron_lengths$count) / sum(intron_lengths$count)
intron_lengths$count = cum_count

intron_lengths <- intron_lengths[intron_lengths$length %% (1024 * 8)  == 0 & intron_lengths$length >= 1024 * 8,]
print(intron_lengths)

my_plot <- ggplot(data = intron_lengths, aes(x = log2(length), y = count)) +
	geom_point() +
	ylab("cumulative")
print(my_plot)

dev.off()

theme_set(previous_theme)
