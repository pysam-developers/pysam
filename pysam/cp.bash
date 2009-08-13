c="bgzf.c kstring.c bam_aux.c bam.c bam_import.c sam.c bam_index.c bam_pileup.c bam_lpileup.c bam_md.c glf.c razf.c faidx.c knetfile.c bam_sort.c  "
h="bam.h bam_endian.h bam_maqcns.h bgzf.h faidx.h glf.h khash.h knetfile.h kseq.h ksort.h kstring.h razf.h sam.h"

for f in ${c} ${h}; do 
    cp ../samtools-0.1.5c/${f} .
done
