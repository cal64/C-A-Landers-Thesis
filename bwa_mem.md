### Designer: Julien Bauer

### Description

Two Python scripts used for alignment of MeDIP-Seq reads

The scripts were called using the following command lines, with inputs and outputs modified as appropriate

```Python
srun python bin/step6_create_reference_indexes.py --ref_path /mnt/cgs-fs3/Sequencing/Genome/Pig/Sus_scrofa.Sscrofa10.2.dna.chromosomes.fa &>reference_indexing_log.txt

srun python bin/step7_mappingReads.py --in_dir trimmedReads/ --out_dir /alignedReads --ref  /mnt/cgs-fs3/Sequencing/Genome/Pig/Sus_scrofa.Sscrofa10.2.dna.chromosomes.fa &>step7_mapping_log.txt
```

### Code

#### Part One: step6_create_reference_indexes.py
Used for indexing of a provided reference genome

```Python

import argparse
import os

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description = 'Prepare the reference for the mapping step by creating indexes')
	parser.add_argument('--ref_path', help='Path to reference file')
	args=parser.parse_args()
	reference = args.ref_path
	reference = reference.strip(".fa")
	os.system("java -jar /usr/local/picard/picard.jar CreateSequenceDictionary REFERENCE="+args.ref_path+" OUTPUT="+reference+".dict")
	os.system("samtools faidx "+args.ref_path)
	os.system("bwa index "+args.ref_path)

```

#### Part Two: step7_mappingReads.py
Used for calling BWA-MEM and completeing alignment

```Python
import os
import sys
import re
import errno
import glob
import time
import pickle
import logging
from joblib import Parallel, delayed
import argparse
import multiprocessing




def make_sure_path_exists(path):
    try:
        os.makedirs(path)
    except OSError as exception:
        if exception.errno != errno.EEXIST:
            raisec

def get_filepaths(directory):
    """
    This function will generate the file names in a directory
    tree by walking the tree either top-down or bottom-up. For each
    directory in the tree rooted at directory top (including top itself),
    it yields a 3-tuple (dirpath, dirnames, filenames).
    """
    file_paths = []  # List which will store all of the full filepaths.

    # Walk the tree.
    for root, directories, files in os.walk(directory):
        for filename in files:
            # Join the two strings in order to form the full filepath.
            filepath = os.path.join(root, filename)
            file_paths.append(filepath)  # Add it to the list.

    return file_paths

def alignment(i):
    if R2exist:
        #count+=1
        r1 = args.in_dir+"/"+i
        r2 = r1.replace("R1","R2")
        r2 = r2.replace("val_1", "val_2")
        mapped=r1.replace("_R1","")
        mapped=mapped.replace(".fq","")
        mapped=mapped.replace(args.in_dir,args.out_dir+"/")
       # print(mapped)
        #print(r1+" "+r2)
        isGZ = re.compile("fq.gz")
        r1IsGZ = isGZ.findall(r1)
        r2IsGZ = isGZ.findall(r2)
        if r1IsGZ is not None:
            print("Unzipping r1 "+r1)
            os.system("gunzip "+r1)
            r1 = r1.strip(".gz")

        if r2IsGZ is not None:
            print("Unzipping r2 "+r2)
            os.system("gunzip "+r2)
            r2 = r2.strip(".gz")

        os.system("bwa mem -t 8 "+ refGenome +" "+r1 + " "+r2 + " >"+mapped+".sam")
       # sys.exit("Mapping done")
        os.system("samtools view -b "+mapped+".sam"+" >"+mapped+".bam")
        os.system("samtools sort -n -o "+mapped+".sortedQuery.bam "+mapped+".bam")
        os.system("samtools sort -o "+mapped+".sortedCoor.bam "+mapped+".bam")
        os.system("samtools index "+mapped+".sortedQuery.bam")
        os.system("samtools index "+mapped+".sortedCoor.bam")

    else:
        print("Error this pipeline needs pair end reads!")




###################
#  * Rscripts used:
# mapping_summary.R
# mapping_distribution.R
###################

parser = argparse.ArgumentParser(description = 'Extract sequences from fastQ file with custom indexes at the beginning of the reads')
parser.add_argument('--in_dir', help='Path to folder containing fastq files. Default=rawReads/', default='rawReads/')
parser.add_argument('--out_dir', help='Path to out put folder. Default=processedReads/', default='processeddReads/')
parser.add_argument('--ref', help='Path to the reference genome')
args=parser.parse_args()

logFilename = './' + sys.argv[0].split(".")[0].split('/')[-1]
logging.basicConfig(format='%(asctime)-5s %(message)-10s ',
                    filename=logFilename + ".log",
                    #filemode='w',
                    level=logging.DEBUG,
                    datefmt='%m-%d-%Y  %H:%M:%S')

logging.info(" ")
logging.info(" ")
logging.info("***************************************")
logging.info("*************MAPPING READS*************")
logging.info(" ")
logging.info("User command: " + str(sys.argv))


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... READING PARAMETERS FILE")
logging.info('-- Input parameters:')
path = args.in_dir
logging.info('Working Directory = ' + path)
refGenome = args.ref
logging.info('Reference Genome = ' + refGenome)



logging.info(" ")
logging.info(" ")
logging.info("#################################")
# Read in sampleNames
logging.info('... Sample names:')
#fastQtoMap = []
fastQFile=os.listdir(args.in_dir)
R1tab=[]
R2exist=False
for i in fastQFile:
    findR1 = re.compile('.+\_R1\_.+.fq')
    findR2 = re.compile('.+\_R2\_.+.fq')
    R1=findR1.match(i)
    R2=findR2.match(i)
    #print(R1)
    if R1 is not None:
        R1tab.append(i)
    if R2 is not None:
        R2exist=True


fastQtoMap= [i for i, x in enumerate(fastQFile) if re.findall("fq", x)]
print(R1tab)
#print(R2exist)


logging.info(" ")
logging.info(" ")
logging.info("#################################")
logging.info("... Mapping of reads")
make_sure_path_exists(args.out_dir)
alignmentOut = os.listdir(args.out_dir) # Check if aligned reads files exist
indicesAlignmentFiles = [i for i, x in enumerate(alignmentOut) if re.findall(".sam", x)] # Check if aligned reads files exist

os.system("mkdir "+args.out_dir+"/sorted")

if not indicesAlignmentFiles:
    Parallel(n_jobs=1)(delayed(alignment)(i) for i in R1tab)
logging.info("... Finished Mapping")
logging.info("Post mapping clean up")
os.system("rm "+args.out_dir+"/*sam")
os.system("mv "+args.out_dir+"/*sorted*bam "+args.out_dir+"/sorted/")
#os.system("rm "+args.out_dir+"/*bam")
#os.system("mv "+args.out_dir+"/sorted/*bam "+args.out_dir)
#os.system("rmdir "+args.out_dir+"/sorted/")
logging.info("Done!")

```
