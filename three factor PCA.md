### Designer: Maria Gomez

### Description

Fragment of code designed to be inserted into the CGS pipeline differential expression script, in order to produce a PCA plot displaying PCs 1 through 3

### Code

```R

PCA
print('Plotting PCA')
data_rlog=assay(rld) # make matrix with rlog values
pca_data=plotPCA(data_rlog, returnData=TRUE, number_of_pcs=number_of_pcs)
percentVar=round(100*attr(pca_data,'percentVar'))
rownames(sample_info)=sample_info$SampleID
pca_data=cbind(pca_data,sample_info[rownames(pca_data),])

pca_data$pch=1
pch=read.table('legend_pca.txt')
pca_data$pch=pch[,1]
pca_data$pch[pca_data$Agression=='Case' & pca_data$Parity==1]=15
pca_data$pch[pca_data$Agression=='Case' & pca_data$Parity==2]=17
pca_data$pch[pca_data$Agression=='Case' & pca_data$Parity==3]=19
pca_data$pch[pca_data$Agression=='Control' & pca_data$Parity==1]=0
pca_data$pch[pca_data$Agression=='Control' & pca_data$Parity==2]=2
pca_data$pch[pca_data$Agression=='Control' & pca_data$Parity==3]=1

pdf(paste0(outdir,'/pca.pdf'), width=11)
pairs(pca_data[,c(1:number_of_pcs)], col=pca_data$Breed, pch=pca_data$pch, oma=c(4,4,4,16), cex=1.5, labels=paste(c(paste0('PC',1:number_of_pcs)), paste0(percentVar,'%'),sep='\n'))
legend('topright',legend=c('HampxLW', 'Case-P1', 'Case-P2', 'Case-P3', 'Control-P1','Control-P2','Control-P3', '', 'LW', 'Case-P1', 'Case-P2', 'Case-P3', 'Control-P1','Control-P2','Control-P3'), col=c(rep('black',8), rep('red',8)), pch=c(NA, 15,17,20,0,2,1,NA,NA,15,17,20,0,2,1), inset=c(0,0), xpd=TRUE)
dev.off()

```
