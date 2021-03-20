# We could try to do a full assembly of all individual libraries but that would take quite a lot of resources and
# might not even be that great. One way to approach this type of analysis is to normalize the whole dataset down
# to a maximum coverage of individual kmers using khmer. We're going to follow large sections of the The Eel Pond mRNAseq Tutorial 
# described here: https://khmer-protocols.readthedocs.io/en/v0.8.2/mrnaseq/index.html. We need to have the khmer scripts in our
# environment

## NOTE: it is very important that your shell environment is setup in a manner described by the above khmer link!

# We've already finished step one of the Eel Pond workflow, let's move to step two where we catalog kmers of the full
# dataset across all the different treatments and replicates

normalize-by-median.py -p -k 20 -C 20 -N 4 -x 3e9 --savehash normC20k20.kh *.pe.fq.gz
normalize-by-median.py -C 20 --loadhash normC20k20.kh --savehash normC20k20.kh *.se.fq.gz

# Next we filter likely erroneous kmers

filter-abund.py -V normC20k20.kh *.keep

# Now we need to do a bit of renaming and also figure out what is still paired and what has been orphaned again

for i in *.pe.*.abundfilt;
do
   extract-paired-reads.py $i
done

# combine orphaned reads

for i in *.se.fq.gz.keep.abundfilt
do
   pe_orphans=$(basename $i .se.fq.gz.keep.abundfilt).pe.fq.gz.keep.abundfilt.se
   newfile=$(basename $i .se.qc.fq.gz.keep.abundfilt).se.keep.abundfilt.fq.gz
   cat $i $pe_orphans | gzip -c > $newfile
done

# a bit of renaming

for i in *.abundfilt.pe
do
   newfile=$(basename $i .fq.gz.keep.abundfilt.pe).keep.abundfilt.fq
   mv $i $newfile
   gzip $newfile
done

# and a bit of cleaning

rm *.se.fq.gz.keep.abundfilt
rm *.pe.fq.gz.keep.abundfilt.se
rm *.keep
rm *.abundfilt

## be careful here with removing files that don't have the qc suffix

# Now we can use Trinity to do a de novo assembly of our full, normalized dataset
## NOTE: Trinity has changed a lot over the years. This is a very generic command but 
# it works!

Trinity.pl --left left.fq --right right.fq --seqType fq -JM 250G

# We have a transcriptome! Let's move to understanding the expression of different clones in our panel and compare
# how these transcripts localize to the PR locus in Iinb1 vs. Xinb3.