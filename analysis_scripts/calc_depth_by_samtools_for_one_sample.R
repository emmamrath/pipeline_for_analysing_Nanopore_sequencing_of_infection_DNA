options(width=250)
options(max.print=10000)
options(stringsAsFactors = FALSE)
options(scipen=10000)
library(data.table)
library(dplyr)
library(tidyr)
library(stringr)
library(ggplot2)
library(sqldf)
library(reshape2)
library(scales)

args = commandArgs(trailingOnly=TRUE) # for production
#args=c( "ERR3077518_2697019to762948.samtools_depth.head10000.txt", "ERR3077518_2697019to762948.contig_names.head10000.txt", "ERR3077518", "ERR3077518_2697019to762948.samtools_depth.summary.txt" ) # for testing
depth_file = as.character(args[1])
contig_names_file = as.character(args[2])
sample_id = as.character(args[3])
output_file = as.character(args[4])

depth = read.table( depth_file, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
colnames(depth) = c("ncbi_id", "position", "depth")
depth$position = as.numeric(as.character(depth$position))
depth$depth = as.numeric(as.character(depth$depth))

contig_names = read.table( contig_names_file, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
colnames(contig_names) = c("ncbi_id", "contig_name", "contig_length")
contig_names$contig_length = as.numeric(as.character(contig_names$contig_length))

names_depths = merge( x=contig_names, y=depth, by=c("ncbi_id"), all.x=TRUE, all.y=TRUE )

names_depths_ge10 = names_depths[(names_depths$depth>=10),]
names_depths_ge100 = names_depths[(names_depths$depth>=100),]
names_depths_ge1000 = names_depths[(names_depths$depth>=1000),]

names_depths_summary = sqldf("select ncbi_id, contig_name, contig_length, avg(depth), count(*) from names_depths group by ncbi_id, contig_name")
names(names_depths_summary)[names(names_depths_summary) == 'avg(depth)'] <- 'avg_depth'
names(names_depths_summary)[names(names_depths_summary) == 'count(*)'] <- 'num_nucleotides'
names_depths_summary$percent_nucleotides = names_depths_summary$num_nucleotides / names_depths_summary$contig_length * 100

names_depths_summary_ge10 = sqldf("select ncbi_id, contig_name, contig_length, avg(depth), count(*) from names_depths_ge10 group by ncbi_id, contig_name")
names(names_depths_summary_ge10)[names(names_depths_summary_ge10) == 'avg(depth)'] <- 'avg_depth_ge10'
names(names_depths_summary_ge10)[names(names_depths_summary_ge10) == 'count(*)'] <- 'num_nucleotides_depth_ge10'
names_depths_summary_ge10$percent_nucleotides_depth_ge10 = names_depths_summary_ge10$num_nucleotides_depth_ge10 / names_depths_summary_ge10$contig_length * 100

names_depths_summary_ge100 = sqldf("select ncbi_id, contig_name, contig_length, avg(depth), count(*) from names_depths_ge100 group by ncbi_id, contig_name")
names(names_depths_summary_ge100)[names(names_depths_summary_ge100) == 'avg(depth)'] <- 'avg_depth_ge100'
names(names_depths_summary_ge100)[names(names_depths_summary_ge100) == 'count(*)'] <- 'num_nucleotides_depth_ge100'
names_depths_summary_ge100$percent_nucleotides_depth_ge100 = names_depths_summary_ge100$num_nucleotides_depth_ge100 / names_depths_summary_ge100$contig_length * 100

names_depths_summary_ge1000 = sqldf("select ncbi_id, contig_name, contig_length, avg(depth), count(*) from names_depths_ge1000 group by ncbi_id, contig_name")
names(names_depths_summary_ge1000)[names(names_depths_summary_ge1000) == 'avg(depth)'] <- 'avg_depth_ge1000'
names(names_depths_summary_ge1000)[names(names_depths_summary_ge1000) == 'count(*)'] <- 'num_nucleotides_depth_ge1000'
names_depths_summary_ge1000$percent_nucleotides_depth_ge1000 = names_depths_summary_ge1000$num_nucleotides_depth_ge1000 / names_depths_summary_ge1000$contig_length * 100

m1 = merge( x=names_depths_summary, y=names_depths_summary_ge10, by=c("ncbi_id","contig_name","contig_length"), all.x=TRUE, all.y=TRUE )
m2 = merge( x=names_depths_summary_ge100, y=names_depths_summary_ge1000, by=c("ncbi_id","contig_name","contig_length"), all.x=TRUE, all.y=TRUE )
data_summary = merge( x=m1, y=m2, by=c("ncbi_id","contig_name","contig_length"), all.x=TRUE, all.y=TRUE )

data_summary[is.na(data_summary)] <- 0
data_summary$sample_id = sample_id
data_summary = data_summary[,c( ncol(data_summary), 1:(ncol(data_summary)-1) )]

write.table( data_summary, file=output_file, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )


