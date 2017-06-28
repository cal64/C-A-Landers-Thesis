### Designer: Maria Gomez

### Description 
Removes all reads mapping to known ribosomal RNA genes from post-counting .txt files

### Code

#### Part One
To retrieve all the lines from the GTF file which are both described as genes and labeled as rRNA in the gene biotype field. 
These lines are printed to a new file, from which everything except for the ensembl codes of the genes is removed
```python
grep 'rRNA' Sus_scrofa.Sscrofa10.2.85.gtf | awk '{ if($3=="gene") print $0}' | awk '{print $10}' | sed 's/"//g' | sed 's/;//g' > rRNA_genes.txt
```

#### Part Two
To remove rRNA lines from post-counting BAM files using the text file created above

```python
mkdir rRNA_removed ###create new directory for cleaned files

for i in $(cat ../sample_names.txt); do grep -vFf /mnt/cgs-fs3/Sequencing/Genome/Pig/ensembl/gtf/Sscrofa10.2/release-85/rRNA_genes.txt ${i}_count.txt > rRNA_removed/${i}_count.txt ; done ###retrieves all the lines that match rRNA reads, in a loop, for all samples

cd rRNA_removed

mkdir check_order ###create new directory to check gene order within count files

for i in $(cat ../../../sample_names.txt); do awk '{print $1}' ${i}_count.txt > check_order/${i}_count.txt; done ###copies count files into new directory to check gene order, in a loop, for all samples

cd check_order

###lines below use R to check the order of the genes has not changed###

# OPEN R
indir='./'
files=list.files(indir, pattern = "_count.txt", full.names=TRUE)

data=read.table(files[1], stringsAsFactors=FALSE)
colnames(data)=files[1]
rownames(data)=data[,1]

for (i in 2:length(files)){
  name=files[i]
  assign(name, read.table(name, stringsAsFactors=FALSE))
  sample=eval(as.name(name))
  #index=match(rownames(data), sample[,1])
  data=data.frame(data, sample[,1])
  data[,i]=as.character(data[,i])
     }

```
