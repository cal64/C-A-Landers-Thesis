### Designers: 

### Description


Arguments file - reference pipeline


### Code

```R

suppressMessages(library(DESeq2))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(matrixStats))
suppressMessages(library(calibrate))


######################################################################
                              Setup
######################################################################

# Read arguments from file #
#--------------------------#

arguments_file ='arguments_files/deseq2_arguments_basic_WSS_aggressionplussibship.txt' 

arguments_file=read.csv(arguments_file, header=FALSE)

read_args_line=function(x){
  x=gsub(' *$', '',x)
  x=gsub('.*=', '',x)
  x=gsub('^ ','',x)
  return(x)
}
indir=read_args_line(grep('^indir =', arguments_file$V1, value=TRUE))
outdir=read_args_line(grep('^outdir =', arguments_file$V1, value=TRUE))
sample_info=read_args_line(grep('^sample_info =', arguments_file$V1, value=TRUE))
comparisons=read_args_line(grep('^comparisons =', arguments_file$V1, value=TRUE))
design=read_args_line(grep('^design =', arguments_file$V1, value=TRUE))
gtf.file=read_args_line(grep('^gtfFile =', arguments_file$V1, value=TRUE))
main_factor=read_args_line(grep('^main_factor =', arguments_file$V1, value=TRUE))
number_of_pcs=as.numeric(read_args_line(grep('^number_of_pcs =', arguments_file$V1, value=TRUE)))


# Reading counted reads files #
#-----------------------------#

files=list.files(indir, pattern = "_count.txt", full.names=TRUE) # will search for files within the indir with the pattern 'count.txt'

# Create a matrix of counted data #
#---------------------------------#

data=read.table(files[1], row.names=1) # read the first sample count file into data frame
colnames(data)=files[1] 

# Reads the other sample count files in a for loop 
for (i in 2:length(files)){ 
  name=files[i]
  assign(name, read.table(name))
  sample=eval(as.name(name))
  index=match(rownames(data), sample[,1]) # get the samples in the same order 
  data=cbind(data, sample[index,2]) # Add sample to data frame
  colnames(data)[i]=name # change colname of current sample
} 

colnames(data)=gsub('_count.txt','',colnames(data)) # change colnames to remove '_count.txt' suffix from filenames
colnames(data)=gsub('.*/','',colnames(data)) # change colnames to remove everything before the '/'
total_counts=colSums(data) # add the number of counts in each column/sample
print(paste0('Matrix of counts created, \nSamples: ', paste(colnames(data), collapse=' '), ' \nTotal Number of counts:', paste(total_counts, collapse=' ')))

## Remove counts not mapped to features ##
data=data[-c(which(rownames(data)=='__no_feature'):dim(data)[1]),] #removing rows from matrix 
total.data.counted=colSums(data) # add the number of counts in each column/sample

# Getting information about the groups #
#--------------------------------------#

print('Reading group info')
print(sample_info)
sample_info = read.csv(sample_info)
column_main_factor=which(colnames(sample_info)==main_factor)
sample_info=sample_info[order(sample_info[,column_main_factor]),] # order samples according to aggression

# Matching the data matrix to the sample_info table #
#---------------------------------------------------#

data=data[,match(as.character(sample_info$SampleID),colnames(data))] # order data frame based on sample_info order

# Creating a DDS object to analyse data with DESEQ2 #
#---------------------------------------------------#
print('Creating a DDS object')
colData<-data.frame(Group=sample_info$Aggression, Breed=sample_info$Breed, Parity=sample_info$Parity, Sibship=sample_info$Sibship) # reads variable info for each sample into data frame
colData$Parity=as.factor(colData$Parity)
rownames(colData)=sample_info$SampleID
comparisons=read.csv(comparisons)
colData$Group=relevel(colData$Group, ref=as.character(comparisons$baselineGroup)) # specify baseline variable
print(colData)

###the following sections ‘actually’ create the dds object. Variations on the model are provided; only one option can be used at a time###

#Model: Group Only
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Group)###factor desired in results must be last
print(dds)

#Model: Breed+Group
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Breed+Group)###factor desired in results must be last
print(dds)

#Model: Parity+Group
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Parity+Group)###factor desired in results must be last
print(dds)

#Model: Parity+Breed+Group
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Parity+Breed+Group)###factor desired in results must be last
print(dds)

#Model: Sibship+Group
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Sibship+Group)###factor desired in results must be last
print(dds)

# Remove genes with zero counts in all samples  #
#-----------------------------------------------#

dds <- dds[rowSums(counts(dds)) >= 1, ]
dds <- DESeq(dds)
print(dds)

print('DDS object created and pre-filtered. Ready to explore the data')

######################################################################
                           Quality Control
######################################################################

# Transforming values to rlog #
#-----------------------------#

print('Transforming values to rlog')
rld <- rlog(dds, blind=FALSE)
#rld_blind <- rlog(dds, blind=TRUE)

# Library Size #
#--------------#

## Define order of samples in dds ##
sample_order <- c("1299_S5", "1301_S3", "1313_S2", "1311_S1", "358_S7", "2167_S8", "365_S4", "2169_S6", "5219_S7", "5221_S3", "5572_S2", "3028_S1", "5585_S8", "5587_S6", "5845_S5", "5844_S4") #creates list with set order of sample names
ddslib <- dds[,sample_order] # sorts the dds using that list
print(ddslib) # prints dds object to confirm sample order has been changed 

## Remove _SX from sample names ##
colnames(ddslib)<- c("1299", "1301", "1313", "1311", "358", "2167", "365", "2169", "5219", "5221", "5572", "3028", "5585", "5587", "5845", "5844") # manually specifies sample names
colnames(ddslib) # prints dds object to confirm sample names changed

## Create Libsize plot ##
pdf(paste0(outdir,'/Libsize', ".pdf"), width=9)
par(mar=c(12.8,7.1,4.1,2.1), mgp=c(5,1,0))
barplot(colSums(counts(ddslib, normalize=FALSE)), names=colnames(counts(ddslib, normalize=FALSE)), las=2, ylab="Counts", xlab="Samples", col = c("#a6cee3","#a6cee3","#1f78b4","#1f78b4","#b2df8a","#b2df8a","#33a02c","#33a02c","#e31a1c","#e31a1c","#fb9a99","#fb9a99","#fdbf6f","#fdbf6f","#ff7f00","#ff7f00")) # uses same colour palette as PCA
dev.off()

# Heatmap #
#---------#

## Define order of samples in rld ##
colnames(rld)
sample_order <- c("1299_S5", "1301_S3", "1313_S2", "1311_S1", "358_S7", "2167_S8", "365_S4", "2169_S6", "5219_S7", "5221_S3", "5572_S2", "3028_S1", "5585_S8", "5587_S6", "5845_S5", "5844_S4") #creates list with set order of sample names
rldheat <- rld[,sample_order] # sorts rld using that list
colnames(rldheat) # prints rld to confirm sample order has been changed 

## Remove _SX from sample names in rldheat ##
colnames(rldheat)<- c("1299", "1301", "1313", "1311", "358", "2167", "365", "2169", "5219", "5221", "5572", "3028", "5585", "5587", "5845", "5844")
colnames(rldheat) # prints dds object to confirm sample names changed 

## Create Heatmap ##
heatmap_palette <- colorRampPalette(brewer.pal(11,"RdYlGn"))(n = 100) # defines colour palette for heatmap
png(paste0(outdir,'/Heatmap',".png"), width = 5*300, height = 5*300, res = 300, pointsize = 8)
heatmap.2(cor(assay(rldheat)),
  main = "Heatmap",        # creates heatmap title
  notecol="black",         # change font color of cell labels to black
  scale=c('none'),
  density.info="density",  # turns on density plot inside color legend
  trace="none",            # turns off trace lines inside the heat map
  margins =c(12,9),        # widens margins around plot
  col=heatmap_palette,     # use color palette defined earlier
  dendrogram="row")        # only draw a row dendrogram
dev.off()

# Principal Component Analysis #
#------------------------------#

## Calculating PCAs ##
data_rlog=assay(rld)
rv <-rowVars(data_rlog) # takes the variance of all the rows.
ntop = 500
select <- order (rv, decreasing = TRUE) [seq_len(min(ntop,length(rv)))] # reorders vector with most variant gene (and its variance value) at the top
pca<- prcomp(t(data_rlog[select, ])) # takes everything from the data_rlog data frame associated with the top 500 variant genes and transposes the rows and columns
percentVar <- 100*(pca$sdev^2/sum(pca$sdev^2)) # uses the standard deviation of each eigenvector to calculate the variation accounted for by each eigenvector 
d <- data.frame(PC1= pca$x[,1], PC2= pca$x[,2], PC3= pca$x[,3], PC4= pca$x[,4], as.data.frame (colData(dds)), name = colnames(data_rlog)) # creates a new data frame called d, which is sixteen rows )sample names) and four columns (the first four PCs)
print(d)

## Define order of samples in dataframe ##
sample_order <- c("1299_S5", "1301_S3", "1313_S2", "1311_S1", "358_S7", "2167_S8", "365_S4", "2169_S6", "5219_S7", "5221_S3", "5572_S2", "3028_S1", "5585_S8", "5587_S6", "5845_S5", "5844_S4") #creates list with set order of sample names
d1 <- d[match(sample_order, d$name),] #creates new data frame with samples in the same order as sample_order
print(d1) # prints data frame to confirm sample order changed

## Remove _SX from sample names ##
d1$name = gsub("_S.", "", d1$name)
print(d1$name) # prints object to confirm sample names changed 
 
## Creating Plot ##

# manually assigns shapes according to case (15) and controls (19)
# manually assigns colours by control/case pair
# colours to use chosen using http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=8

png(paste0(outdir,'/PCA_plot',".png"), width = 6*300, height = 5*300, res = 300, pointsize = 8) # png version
pairs(d1[,(1:3)], 
  col = c("#a6cee3","#a6cee3","#1f78b4","#1f78b4","#b2df8a","#b2df8a","#33a02c","#33a02c","#e31a1c","#e31a1c","#fb9a99","#fb9a99","#fdbf6f","#fdbf6f","#ff7f00","#ff7f00"), # sets colours for each sample
  oma=c(4,4,4,16), 
  pch= c(19,15,19,15,19,15,19,15,19,15,19,15,19,15,19,15), # sets shape for each sample
  cex=2, 
  labels=paste(c('PC1','PC2','PC3'), # labels which rows and columns of grid correspond to each PC
  paste0(round(percentVar, digits = 2),'%'),sep='\n')) # prints percent variance each PC accounts for
legend('right', legend= d1$name, 
  cex=1.3, 
  col= c("#a6cee3","#a6cee3","#1f78b4","#1f78b4","#b2df8a","#b2df8a","#33a02c","#33a02c","#e31a1c","#e31a1c","#fb9a99","#fb9a99","#fdbf6f","#fdbf6f","#ff7f00","#ff7f00"), 
  pch= c(19,15,19,15,19,15,19,15,19,15,19,15,19,15,19,15), 
  inset=c(0,0), 
  xpd=TRUE)
dev.off()

pdf(paste0(outdir,'/PCA_plot',".pdf"), width=11, height=9) # pdf version 
pairs(d1[,(1:3)], col = c("#a6cee3","#a6cee3","#1f78b4","#1f78b4","#b2df8a","#b2df8a","#33a02c","#33a02c","#e31a1c","#e31a1c","#fb9a99","#fb9a99","#fdbf6f","#fdbf6f","#ff7f00","#ff7f00"), oma=c(4,4,4,16), pch= c(19,15,19,15,19,15,19,15,19,15,19,15,19,15,19,15), cex=2, labels=paste(c('PC1','PC2','PC3'), paste0(round(percentVar, digits = 2),'%'),sep='\n'))
legend('right', legend= d1$name, cex=1.3, col= c("#a6cee3","#a6cee3","#1f78b4","#1f78b4","#b2df8a","#b2df8a","#33a02c","#33a02c","#e31a1c","#e31a1c","#fb9a99","#fb9a99","#fdbf6f","#fdbf6f","#ff7f00","#ff7f00")
, pch= c(19,15,19,15,19,15,19,15,19,15,19,15,19,15,19,15), inset=c(0,0), xpd=TRUE)
dev.off()


######################################################################
                      Differential Expression
######################################################################

# Size factor correlation with total counts #
#-------------------------------------------#

pdf(paste0(outdir,'/Size_factor_vs_counts',".pdf"))
plot(sizeFactors(dds),colSums(counts(dds)))
abline(lm(colSums(counts(dds)) ~ sizeFactors(dds) + 0))
dev.off()


# Transforming values to rlog #
#-----------------------------#

rld <- rlog(dds, blind=FALSE)
    
# Dispersion Plot #
#------------#

png(paste0(outdir,'/ Dispersion',".png"))
plotDispEsts(dds)
dev.off()

# Differential Expression Calculations #
#--------------------------------------#

res <- results(dds)
res <- res[order(res$padj),]
head(res)	
write.csv(as.data.frame(res), paste0(outdir,'/results_data_frame', ".csv"))

# MA Plot #
#---------#

png(paste0(outdir,'/MAplot',".png"))
plotMA(res)
dev.off()

# Load Annotation from GTF #
#--------------------------#

GTF <- import.gff(gtf.file, format="gtf", feature.type="gene") # read gtf file
df=GTF$gene_name # get gene_names
names(df)=GTF$gene_id # match with gene ids
 
no_name=df[-which(df!='NA')] # use ensembl symbol for empty names 
no_name=names(no_name) # use ensembl symbol for empty names
names(no_name)=no_name # use ensembl symbol for empty names
df=df[which(df!='NA')] # get genes with names
df=c(df, no_name) # paste genes with names and genes with empty names(now with ensembl ids)

# Annotate the data frame #
#-------------------------#

res.df <- as.data.frame(res) # creating a data frame
detags <- rownames(res.df) 
res.df = cbind(res.df,counts(dds,normalized=TRUE)[detags,]) # add normalised count values to the data frame and sorts all the genes in the same order as the results data frame
res.df = cbind(GeneName=df[rownames(res.df)], res.df) # adding the gene names to the data frame
write.csv(res.df, paste0(outdir,'/annotated_results',".csv")) # writes data frame to csv file

# Pvalue distribution #
#---------------------#

pvalue_palette <- colorRampPalette(brewer.pal(08,"RdGy"))(n = 200) # specifies colour palette

png(paste0(outdir,'/Pvalue_distribution',".png"), width = 6*300, height = 5*300, res= 300, pointsize = 8) # png version
hist(res$pvalue, xlab='p-value', main='', breaks=100, col=pvalue_palette)
dev.off() 

pdf(paste0(outdir,'/Pvalue_distribution',".pdf"), width = 10, height = 9) # pdf version
hist(res$pvalue, xlab='p-value', main='', breaks=100, col=pvalue_palette)
dev.off() 

# Volcano Plot #
#--------------#

# modified from http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html 

png(paste0(outdir, '/Volcano_plot',".png"))
with(res.df, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
#colored points: red if padj<0.1, orange if <0.05, green if <0.01)
with(subset(res.df, padj<0.1 ), points(log2FoldChange, -log10(pvalue), pch = 20, col="red"))
with(subset(res.df, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch = 20, col="orange"))
with(subset(res.df, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch = 20, col="green"))
# Label points where padj<0.1 with the textxy function from the calibrate plot
with(subset(res.df, padj<0.1), textxy(log2FoldChange, -log10(pvalue), labs=GeneName, cex=.8))
dev.off()

```
