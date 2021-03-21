# We're going to do our DE analysis using the kallisto mapper. Let's start by indexing the de novo transcriptome from Trinity

kallisto index -i transcripts.idx Trinity.fasta

# Now, let's (pseudo-)map all of our different treatments and replicates to the reference
# the basic form is:

kallisto quant -i transcripts.idx -o output -b 100 <(zcat reads_1.fastq.gz) <(zcat reads_2.fastq.gz)

# We're going to use parallel to speed this up. The script has the form:

#!/bin/sh

sample=$1
describer=$(echo ${sample} | sed 's/_R1.fastq.gz//')

kallisto quant -i transcripts.idx -o ${describer} -b 100 <(zcat ${describer}_R1.fastq.gz) <(zcat ${describer}_R2.fastq.gz)

# we can then execute this script with the following:

ls *_R1.fastq.gz | parallel -j24 -k bash script.sh {}

# where script is from lines 12-18
# okay, so we have the counts for each treatment/replicate in a directory. Now we just need to compile the results so that we can use
# sleuth so that we can do the DE analysis. So now we just need to open an R shell and call sleuth of the directory we created. See 
# sleuth.R