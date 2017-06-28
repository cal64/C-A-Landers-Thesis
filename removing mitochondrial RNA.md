### Designer: Maria Gomez

### Description 
Removes all reads mapping to mitochondrial chromosome from post-alignment, sorted BAM files

### Code
```python
samtools idxstats alignedBAMfile.bam | cut -f 1 | grep -v MT | xargs samtools view -b alignedBAMfile.bam > cleanedfiledirectory/alignedandcleanedBAMfile.bam
```
