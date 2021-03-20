# The first thing we need to do is pull the data from stressflea. The full dataset
# can be obtained using either NCBI or EMBL. All the necessary accessions are provided
# in the data directory. Once we have the data we will proceed with quality assessment 
# and filtering

# First, let's run fastqc in the directory which contains all left and right reads for 
# each treatment or replicate

fastqc -t 24 *.fastq.gz

# As expected we've got primer contamination and right read degredation. So let's use
# trimmomatic to remove adapters and do a bit of light trimming. We're going to use a
# parallels hack so that we can run this over the whole directory of raw data. The 
# script needs to be called with the follow command:

ls *_R1.fastq.gz | parallel -j4 -k bash trimmomatic.sh {}

# so now we've trimmed the data, and interleaved paired reads and concatenated orphaned
# reads. Let's repeat the fastqc run and then combine data with multiqc (after removing
# previous fastqc runs from the directory

rm *.zip *.html
fastqc -t 24 *.fq.gz
multiqc *

# things look good! So now we can move towards the de novo assembly of the full transcriptome dataset!





