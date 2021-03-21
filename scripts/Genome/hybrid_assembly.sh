# Below are the necessary commands needed to generate a de novo assembly of the Iinb1 clone. Note, other assembly methods were attempted but 
# these did not result in a complete assembly of the ABC, or PR, locus. These include masurca for short-reads and then cabog/canu for the 
# PacBio data. We have relatively low coverage for the PacBio so it's not surprising that the hybrid assembly was better.

# We're going to use DBG2OLC in order to conduct the assembly. We've already pulled the Illumina data from the DBSSE and PacBio data from 
# the Functional Genomics Center (FGC).

# let's check the contents of the Illumina data

fastqc -t 2 Iinb1_R1.fastq.gz Iinb1_R2.fastq.gz

# Reads look fine. Let's trim off adapters

java -jar trimmomatic-0.35.jar PE Iinb1_R1.fastq.gz Iinb1_R2.fastq.gz Iinb1_forward_paired.fq.gz \
Iinb1_forward_unpaired.fq.gz Iinb1_reverse_paired.fq.gz Iinb1_reverse_unpaired.fq.gz \
ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# we're going to make a few new files to be used for different steps later on. Mainly let's interleave the paired reads
# let's concatenate the orphaned reads, and let's concatenate both of these files to be used for a high confidence contig set
# assembly

seqtk mergepe Iinb1_forward_paired.fq.gz Iinb1_reverse_paired.fq.gz | gzip -3 > Iinb1.pe.fq.gz
cat Iinb1_forward_unpaired.fq.gz Iinb1_reverse_unpaired.fq.gz > Iinb1.se.fq.gz 
cat Iinb1.pe.fq.gz Iinb1.se.fq.gz > Iinb1.qc.fq.gz

# First step is creating a set of high-quality contigs. We'll use SparseAssembler as recommended by the developer of DBG2OLC; note that 
# other de brujn assemblers were tried but these did not result in improved assembly contiguity

gunzip -c Iinb1.qc.fq.gz > Iinb1.qc.fq
SparseAssembler LD 0 k 51 g 15 NodeCovTh 1 EdgeCovTh 0 GS 23000000 f Iinb1.qc.fq

## multiple k values were tried as well. In the end the important file here is titled Contigs.txt

DBG2OLC k 17 AdaptiveTh 0.0001 KmerCovTh 2 MinOverlap 20 RemoveChimera 1 Contigs Contigs.txt f Iinb1.pacbio.fasta

# The target output here is a file called backbone_raw.fasta. Now we need to polish/generate a consensus file

cat Contigs.txt Iinb1.pacbio.fasta > ctg_pb.fasta
sh split_and_run_sparc.sh backbone_raw.fasta DBG2OLC_Consensus_info.txt ctg_pb.fasta ./consensus_dir 2> cns_log.txt

# Now we have Iinb1_dbg2olc.fasta. Let's try to polish things up a bit more with quiver and then pilon
# First, we map the PacBio data back to our new reference with blasr

blasr Iinb1.subreadset.xml Iinb1_dbg2olc.fasta --bam --out alignments.bam

# Now, let's use quiver to get a polished reference

quiver alignment.bam -r Iinb1_dbg2olc.fasta -o Iinb1_dbg2olc.quiver.fasta

# We've got pretty low coverage PacBio (>15X) so the polish here probably isn't great. Let's add the Illumina
# data in with pilon. First, map Illumina data

bwa index Iinb1_dbg2olc.quiver.fasta
bwa mem -t 55 -p Iinb1_dbg2olc.quiver.fasta Iinb1.qc.fq.gz | samtools sort -O BAM -o Iinb1_dbg2olc.quiver.align.bam -
samtools index Iinb1_dbg2olc.quiver.align.bam

# now pilon

java -jar pilon.jar --genome Iinb1_dbg2olc.quiver.fasta --frags Iinb1_dbg2olc.quiver.align.bam 
mv pilon.fasta Iinb1_dbg2olc.polish.fasta

# We have our new assembly! Now let's find the PR-locus and compare it to the Xinb3 PR-locus!