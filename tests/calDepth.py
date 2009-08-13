import samtools

samtools.test_pileup( "ex1.bam" )

samtools.test_fetch(  "ex1.bam", "seq1:10-20" )
