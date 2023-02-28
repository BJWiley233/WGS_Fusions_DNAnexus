set -o pipefail
set -o errexit

cram=$1
reference=$2
fusion_sites=$3
eid_nameroot=$4

#link the full cram file - necessary because crams and crais can have different bucket addresses
# ln -s ~{cram} allreads.cram
# ln -s ~{cram_index} allreads.crai
# ln -s ~{cram_index} allreads.cram.crai
#pull metrics
samtools flagstat $cram > "$eid_nameroot.allreads.flagstat"
#subset to just the regions of interest
samtools view -T $reference -H $cram -@8 > tmp.sam
cat $fusion_sites | while read chr start stop ; do samtools view -T $reference $cram $chr:$start-$stop >> tmp_reads.sam ; done
sort -S 8G tmp_reads.sam | uniq >> tmp.sam
samtools sort -@8 -O bam -o "$eid_nameroot.subset.bam" tmp.sam
samtools index "$eid_nameroot.subset.bam"
#sanity check for read lengths
samtools view -@8 "$eid_nameroot.subset.bam" | cut -f 10 | perl -nae 'print length($_) . "\n"' | sort | uniq -c | sort -nrk 1 > "$eid_nameroot.read_lengths.txt"
#filter to just high-quality reads
perl /opt/filter_bam_reads.pl "$eid_nameroot.subset.bam" | samtools view -@8 -Sb -o "$eid_nameroot.subset_filtered.bam" -
samtools index "$eid_nameroot.subset_filtered.bam"
rm tmp.sam
rm tmp_reads.sam

