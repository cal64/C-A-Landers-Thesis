### Designers: C A Landers and Maria Gomez

### Description

Code used to perform differential methylation analysis, and produce plots discussed in Chapter Four.
Note that this is not a script - each element must be copied into R and run 'manually'.
For detail on how arguments file is created see https://github.com/CGSbioinfo/RNASeq_pipeline .

### Code

```R

suppressMessages(library(DESeq2))
suppressMessages(library(gplots))
suppressMessages(library(rtracklayer))
suppressMessages(library(ggplot2))
suppressMessages(library(RColorBrewer))
suppressMessages(library(matrixStats))
suppressMessages(library(calibrate))
suppressMessages(library(dplyr))

######################################################################
                              Setup
######################################################################

# Read arguments from file #
#--------------------------#

arguments_file ='arguments_files/arguments_bothtissues_WSS_aggression_breed_parity_tissue.txt' 

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
#main_factor=read_args_line(grep('^main_factor =', arguments_file$V1, value=TRUE))
#number_of_pcs=as.numeric(read_args_line(grep('^number_of_pcs =', arguments_file$V1, value=TRUE)))


# Reading counted reads files #
#-----------------------------#

files=list.files(indir, pattern = "_count.txt", full.names=TRUE)  # will search for files within the indir with the pattern 'count.txt'

# Create a matrix of counted data #
#---------------------------------#

data=read.table(files[1], row.names=1) # read the first sample count file
colnames(data)=files[1]  

# Read the other sample count files in a for loop 
for (i in 2:length(files)){
  name=files[i]
  assign(name, read.table(name))
  sample=eval(as.name(name))
  index=match(rownames(data), sample[,1]) # get the samples in the same order 
  data=cbind(data, sample[index,2]) # Add sample to data (previously created)
  colnames(data)[i]=name # change colname of current sample
}
colnames(data)=gsub('_count.txt','',colnames(data)) # change colnames - remove '_count.txt'
colnames(data)=gsub('.*/','',colnames(data)) # change colnames - remove everything before the '/'
total_counts=colSums(data) # add the number of counts in each column/sample
print(paste0('Matrix of counts created, \nSamples: ', paste(colnames(data), collapse=' '), ' \nTotal Number of counts:', paste(total_counts, collapse=' ')))

print('Keep counts falling to features:')
data=data[-c(which(rownames(data)=='__no_feature'):dim(data)[1]),] #removing rows from data which describe counts that have not been mapped to features
total.data.counted=colSums(data)

# Getting information about the groups #
#--------------------------------------#

print('Reading group info')
print(sample_info)
sample_info = read.csv(sample_info)

# Matching the data matrix to the sample_info table #
#---------------------------------------------------#

data=data[,match(as.character(sample_info$SampleID),colnames(data))] # order data frame based on sample_info order

# Creating a DDS object to analyse data with DESEQ2 #
#---------------------------------------------------#
print('Creating a DDS object')
colData<-data.frame(Group=sample_info$Aggression, Breed=sample_info$Breed, Tissue=sample_info$Tissue, Parity=sample_info$Parity, Sibship=sample_info$Sibship)
colData$Parity=as.factor(colData$Parity)
#colData<-data.frame(Group=sample_info[,column_main_factor])
rownames(colData)=sample_info$SampleID
comparisons=read.csv(comparisons)
colData$Group=relevel(colData$Group, ref=as.character(comparisons$baselineGroup)) # specify baseline variable
print(colData)

#the following sections ‘actually’ create the dds object. Variations on the model are provided; only one option can be used at a time#

#Model: Aggression Only
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Group)###factor desired in results must be last
print(dds)

#Model: Breed+Aggression
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Breed+Group)###factor desired in results must be last
print(dds)

#Model: Parity+Aggression
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Parity+Group)###factor desired in results must be last
print(dds)

#Model: Tissue+Aggression
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Tissue+Group)###factor desired in results must be last
print(dds)

#Model: Parity+Breed+Aggression
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Parity+Breed+Group)###factor desired in results must be last
print(dds)

#Model: Tissue+Breed+Aggression
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Tissue+Breed+Group)###factor desired in results must be last
print(dds)

#Model: Tissue+Parity+Aggression
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Tissue+Parity+Group)###factor desired in results must be last
print(dds)

#Model: Tissue+Parity+Breed+Aggression
dds <- DESeqDataSetFromMatrix(countData= data, colData=colData, design=~Tissue+Parity+Breed+Group)###factor desired in results must be last
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
sample_order <- c("1299H-_S5", "1301H-_S1", "1313H-_S2", "2231H-_S4", "358H-_S3", "358C-Cortex-MethylCollector_S5", "365C-Cortex-MethylCollector_S2", "2169C-Cortex-MethylCollector_S1", "5221C-Cortex-MethylCollector_S4", "4557C-Cortex-MethylCollector_S3")#creates list with set order of sample names
ddslib <- dds[,sample_order]
print(ddslib) 


## Remove _SX from sample names ##
colnames(ddslib)
colnames(ddslib)<- c("1299H", "1301H", "1313H", "2231H", "358H", "358C", "365C", "2169C", "5221C", "4557C")
colnames(ddslib) # prints dds object to confirm sample names changed

## Create Libsize plot ##
pdf(paste0(outdir,'/Libsize', ".pdf"), width=9)
par(mar=c(12.8,7.1,4.1,2.1), mgp=c(5,1,0))
barplot(colSums(counts(ddslib, normalize=FALSE)), names=colnames(counts(ddslib, normalize=FALSE)), las=2, ylab="Counts", xlab="Samples", col = c("#a6cee3","#a6cee3","#1f78b4","#cab2d6","#b2df8a","#b2df8a","#33a02c","#33a02c","#e31a1c","#6a3d9a")) # uses same colour palette as PCA
dev.off()

# Heatmap #
#---------#

## Define order of samples in rld ##
colnames(rld)
sample_order <- c("1299H-_S5", "1301H-_S1", "1313H-_S2", "2231H-_S4", "358H-_S3", "358C-Cortex-MethylCollector_S5", "365C-Cortex-MethylCollector_S2", "2169C-Cortex-MethylCollector_S1", "5221C-Cortex-MethylCollector_S4", "4557C-Cortex-MethylCollector_S3")#creates list with set order of sample names
rldheat <- rld[,sample_order] #use this for Heatmap
colnames(rldheat) # prints rld to confirm sample order has been changed 

## Remove _SX from sample names in rldheat ##
colnames(rldheat)<- c("1299H", "1301H", "1313H", "2231H", "358H", "358C", "365C", "2169C", "5221C", "4557C")
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

## Define order of samples in dataframe##
sample_order <- c("1299H-_S5", "1301H-_S1", "1313H-_S2", "2231H-_S4", "358H-_S3", "358C-Cortex-MethylCollector_S5", "365C-Cortex-MethylCollector_S2", "2169C-Cortex-MethylCollector_S1", "5221C-Cortex-MethylCollector_S4", "4557C-Cortex-MethylCollector_S3")#creates list with set order of sample names
d1 <- d[match(sample_order, d$name),] #creates new data frame with samples in the same order as sample_order
print(d1) 

##remove _SX from sample names###
d1$name = gsub("-.*", "", d1$name)
print(d1$name) #prints dds object to confirm sample names changed 

## Creating Plot ##

# manually assigns shapes according to case (15) and controls (19)
# manually assigns colours by control/case pair
# colours to use chosen using http://colorbrewer2.org/#type=qualitative&scheme=Paired&n=8

png(paste0(outdir,'/PCA_plot',".png"), width = 6*300, height = 5*300, res = 300, pointsize = 8)
pairs(d1[,(1:3)], 
  col = c("#a6cee3","#a6cee3","#1f78b4","#cab2d6","#b2df8a","#b2df8a","#33a02c","#33a02c","#e31a1c","#6a3d9a"), # sets colours for each sample
  oma=c(4,4,4,16), 
  pch= c(19,15,19,15,19,19,19,15,15,15), # sets shape for each sample
  cex=2, 
  labels=paste(c('PC1','PC2','PC3'),  # labels which rows and columns of grid correspond to each PC
  paste0(round(percentVar, digits = 2),'%'),sep='\n')) # prints percent variance each PC accounts for
legend('right', legend= d1$name, 
  cex=1.3, 
  col= c("#a6cee3","#a6cee3","#1f78b4","#cab2d6","#b2df8a","#b2df8a","#33a02c","#33a02c","#e31a1c","#6a3d9a"), 
  pch= c(19,15,19,15,19,19,19,15,15,15), 
  inset=c(0,0), 
  xpd=TRUE)
dev.off()

######################################################################
                      Differential Methylation
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

res=data.frame(unique_id=rownames(res), res) # If res doesnt have “unique_id” column, we add it with this line
GTF <- import.gff(gtf.file, format="gtf", feature.type="exon") # read gtf file
df=data.frame(GTF)
res=left_join(res, df, by='unique_id') # Join res with GTF
write.csv(as.data.frame(res), paste0(outdir,'/results_data_frame', ".csv"))

# Pvalue distribution #
#---------------------#

pvalue_palette <- colorRampPalette(brewer.pal(08,"RdGy"))(n = 200) # specifies colour palette

png(paste0(outdir,'/Pvalue_distribution',".png"), width = 6*300, height = 5*300, res= 300, pointsize = 8) # png version
hist(res$pvalue, xlab='p-value', main='', breaks=100, col=pvalue_palette)
dev.off() 

# Volcano Plot #
#--------------#

# modified from http://www.gettinggeneticsdone.com/2014/05/r-volcano-plots-to-visualize-rnaseq-microarray.html 

png(paste0(outdir, '/Volcano_plot',".png"))
with(res, plot(log2FoldChange, -log10(pvalue), pch=20, main="Volcano plot", xlim=c(-2.5,2)))
#colored points: red if padj<0.1, orange if <0.05, green if <0.01)
with(subset(res, padj<0.1 ), points(log2FoldChange, -log10(pvalue), pch = 20, col="red"))
with(subset(res, padj<.05 ), points(log2FoldChange, -log10(pvalue), pch = 20, col="orange"))
with(subset(res, padj<.01 ), points(log2FoldChange, -log10(pvalue), pch = 20, col="green"))
# Label points where padj<0.1 with the textxy function from the calibrate plot
#with(subset(res, padj<0.1), textxy(log2FoldChange, -log10(pvalue), labs=GeneName, cex=.8))
dev.off()

```
