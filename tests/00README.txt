File ex1.fa contains two sequences cut from the human genome
build36. They were extracted with command:

  samtools faidx human_b36.fa 2:2043966-2045540 20:67967-69550

Sequence names were changed manually for simplicity. File ex1.sam.gz
contains MAQ alignments extracted with:

  (samtools view NA18507_maq.bam 2:2044001-2045500;
   samtools view NA18507_maq.bam 20:68001-69500)

and processed with `samtools fixmate' to make it self-consistent as a
standalone alignment.

In order for the script pysam_test.py to work, you will need pysam
in your PYTHONPATH.

In order for the script example.py to work, you will need pysam
in your PYTHONPATH and run

  make all

beforehand.
