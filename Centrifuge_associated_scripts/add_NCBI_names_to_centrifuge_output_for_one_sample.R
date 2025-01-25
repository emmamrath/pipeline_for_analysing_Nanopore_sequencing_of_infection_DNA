
options(width=200)
options(stringsAsFactors = FALSE)
library(dplyr)
library(data.table)
library(sqldf)

args = commandArgs(trailingOnly=TRUE) # for production
#args=c( 'infile.tsv', 'infile2.tsv', 'outfile.tsv' ) # for testing

infile = as.character(args[1])
inlist = as.character(args[2])
outfile = as.character(args[3])

data = read.table( infile, sep="\t", header=TRUE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
datalist = read.table( inlist, sep="\t", header=FALSE, quote='"', comment.char="", row.names=NULL, stringsAsFactors=FALSE, colClasses=c("character") )
colnames(datalist) = c("seqID", "organism_name")

outdata = merge( x=data, y=datalist, by=c("seqID"), all.x=TRUE, all.y=FALSE )

outdata2 = outdata[ order(outdata$readID, outdata$seqID), c(2,1,3,4,5,6,7,8,9)]

write.table( outdata2, file=outfile, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )

summary1 = sqldf("select seqID, taxID, organism_name, count(*) from outdata2 group by seqID, taxID, organism_name")

names(summary1)[names(summary1) == 'count(*)'] <- 'num_reads'

taxid_summary = sqldf("select taxID, sum(num_reads) from summary1 group by taxID")

names(taxid_summary)[names(taxid_summary) == 'sum(num_reads)'] <- 'num_reads_in_taxID'

summary2 = merge( x=summary1, y=taxid_summary, by=c("taxID"), all.X=TRUE, all.Y=TRUE )

summary2 = summary2[ (summary2$organism_name != "."), ]

summary2 = summary2[ order(summary2$num_reads_in_taxID,summary2$num_reads,decreasing=TRUE), ]

write.table( summary2, file=outfile_summary, quote=FALSE, sep="\t", row.names=FALSE, col.names=TRUE )

