# plot_the_coverage_for_all_contigs_in_a_sample.R

system.file(package="ggplot2")
options(width=250)
library("ggplot2") # for the plot
library("ggrepel") # for spreading text labels on the plot
library("scales") # for axis labels notation

args = commandArgs(trailingOnly=TRUE) # for production
#args=c( 'charalampous2019 ERR3077518', 'align_charalampous2019_reads_to_centrifuge_organisms/sample_genomes/ERR3077518_2697019to762948/ERR3077518_2697019to762948.contig_names.txt', 'align_charalampous2019_reads_to_centrifuge_organisms/minimap2_bam/ERR3077518_2697019to762948.samtools_depth_intervals.head10000.txt', 'align_charalampous2019_reads_to_centrifuge_organisms/minimap2_bam/ERR3077518_2697019to762948.samtools_depth_intervals.head10000.png' ) # for testing

sample = as.character(args[1])
contig_names_file = as.character(args[2])
depth_intervals_file = as.character(args[3])
outfile = as.character(args[4])

outfile_small = paste( gsub(".png$", "", outfile), "_small.png", sep="" )

contig_names = read.table( contig_names_file, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
colnames(contig_names)[1] = 'contig'
colnames(contig_names)[2] = 'contig_name'
colnames(contig_names)[3] = 'size'
contig_names$size = as.numeric(as.character(contig_names$size))

depth_intervals = read.table( depth_intervals_file, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
colnames(depth_intervals)[1] = 'contig'
colnames(depth_intervals)[2] = 'start_pos'
colnames(depth_intervals)[3] = 'end_pos'
colnames(depth_intervals)[4] = 'depth'
depth_intervals$start_pos = as.numeric(as.character(depth_intervals$start_pos))
depth_intervals$end_pos = as.numeric(as.character(depth_intervals$end_pos))
depth_intervals$depth = as.numeric(as.character(depth_intervals$depth))
depth_intervals$feature = "depth"

contig_names = contig_names[(contig_names$contig %in% unique(depth_intervals$contig)),]

contig_sizes = contig_names[,c("contig", "size")]

contig_order = unique(contig_names$contig)
contig_key = setNames(object = as.character(1:length(contig_order)), nm = contig_order)
contig_order = factor(x = contig_order, levels = contig_order)

# convert the contig column in each dataset to the ordered factor
contig_sizes[["contig"]] = factor(x = contig_sizes[["contig"]], levels = contig_order)
depth_intervals[["contig"]] = factor(x = depth_intervals[["contig"]], levels = contig_order)

contig_width_on_plot = 0.1
depth_width_on_plot = 0.6

max_depth = max(depth_intervals$depth)
depth_intervals$depth_size_for_plot = depth_width_on_plot * (depth_intervals$depth / max_depth)
depth_intervals$add_0 = 0
depth_intervals$add_0p1 = ifelse( (depth_intervals$depth_size_for_plot >= 0.1), 0.1, 0)
depth_intervals$add_0p2 = ifelse( (depth_intervals$depth_size_for_plot >= 0.2), 0.2, 0)
depth_intervals$add_0p3 = ifelse( (depth_intervals$depth_size_for_plot >= 0.3), 0.3, 0)
depth_intervals$add_0p4 = ifelse( (depth_intervals$depth_size_for_plot >= 0.4), 0.4, 0)
depth_intervals$add_0p5 = ifelse( (depth_intervals$depth_size_for_plot >= 0.5), 0.5, 0)
depth_intervals$add_0p6 = ifelse( (depth_intervals$depth_size_for_plot >= 0.6), 0.6, 0)

myplot = ggplot(data = contig_sizes) + 
    # base rectangles for the chroms, with numeric value for each chrom on the x-axis
    geom_rect(aes(xmin = as.numeric(contig) - 0, 
                  xmax = as.numeric(contig) - contig_width_on_plot, 
                  ymax = size, ymin = 0),
              colour="black", fill = "white") +
    # black & white color theme 
    theme(axis.text.x = element_text(colour = "black"), 
          axis.text.y = element_text(hjust = 0),
          panel.grid.major = element_blank(), 
          panel.grid.minor = element_blank(), 
          panel.background = element_blank()) + 
    # give the appearance of a discrete axis with chrom labels
    # and remove space between x-axis (will be flipped to become y-axis) and the plot
    scale_x_discrete(name = paste( "contigs, max.depth is ", as.character(max_depth), "X\n(dots show very small intervals)", sep="" ),
          limits = names(contig_key), expand = c(0, 0)) +
    ylab("region on contig (nucleotide position)")

myplot = myplot + 
    # rotate the plot 90 degrees
    coord_flip() +
    ggtitle(paste( sample, " - sequencing depths at contig nucleotide positions", sep="" )) +
    # supress scientific notation on the y-axis
    scale_y_continuous(labels = comma)

myplot = myplot + 
    # add bands for feature value
    geom_rect(data = depth_intervals, aes(xmin = as.numeric(contig) + depth_size_for_plot, 
                                     xmax = as.numeric(contig) + 0, 
                                     ymax = end_pos, ymin = start_pos, fill = feature))

myplot = myplot + 
    # make sure that a pixel point is plotted, even if the interval is too small for a pixel width to be plotted
    geom_point(data = depth_intervals, aes(shape = ".", x = as.numeric(contig) + add_0, y = start_pos, color = feature)) +
    geom_point(data = depth_intervals, aes(shape = ".", x = as.numeric(contig) + add_0p1, y = start_pos, color = feature)) +
    geom_point(data = depth_intervals, aes(shape = ".", x = as.numeric(contig) + add_0p2, y = start_pos, color = feature)) +
    geom_point(data = depth_intervals, aes(shape = ".", x = as.numeric(contig) + add_0p3, y = start_pos, color = feature)) +
    geom_point(data = depth_intervals, aes(shape = ".", x = as.numeric(contig) + add_0p4, y = start_pos, color = feature)) +
    geom_point(data = depth_intervals, aes(shape = ".", x = as.numeric(contig) + add_0p5, y = start_pos, color = feature)) +
    geom_point(data = depth_intervals, aes(shape = ".", x = as.numeric(contig) + add_0p6, y = start_pos, color = feature))

myplot = myplot + theme(legend.position="none")

myplot

# make a smaller summary version
num_contigs = nrow(contig_names)
max_contig_size = max(contig_names$size)
yaxis_width = 4
ggplot_width = yaxis_width + (max_contig_size/1000000 * 1)
header_trailer_height = 1.5
ggplot_height = 0.2 * num_contigs

ggsave(outfile_small, width=ggplot_width, height=ggplot_height, units="cm", limitsize=FALSE)

# adjust plot size according to largest contig and number of contigs, in cm
num_contigs = nrow(contig_names)
max_contig_size = max(contig_names$size)
yaxis_width = 4
ggplot_width = yaxis_width + (max_contig_size/1000000 * 6)
header_trailer_height = 1.5
ggplot_height = 1.5 * num_contigs

ggsave(outfile, width=ggplot_width, height=ggplot_height, units="cm", limitsize=FALSE)

