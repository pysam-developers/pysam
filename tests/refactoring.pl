while (<STDIN>) {
    # Samfile refactoring
    s/Samfile/AlignmentFile/g;

    # AlignedRead refactoring
    s/AlignedRead/AlignedSegment/g;
    
    # Do these patterns first as they match
    # the new names
    s/\.query/\.query_alignment_sequence/g;
    s/\.positions/\.getReferencePositions()/g;

    # basic attributes
    s/\.qname/\.query_name/g;
    s/\.tid/\.reference_id/g;
    s/\.pos/\.reference_start/g;
    s/\.mapq/\.mapping_quality/g;
    s/\.rnext/\.next_reference_id/g;
    s/\.pnext/\.next_reference_start/g;
    s/\.tlen/\.query_length/g;
    s/\.seq/\.query_sequence/g;
    if (/\.qual =/) {
	s/([[\].0-9a-zA-Z]*)\.qual = (\S*)/$1.query_qualities = pysam.fromQualityString($2)/g;
    } else {
	s/([[\].0-9a-zA-Z]*)\.qual/pysam.toQualityString($1\.query_qualities)/g;
    }
    s/\.alen/\.reference_length/g;
    s/\.aend/\.reference_end/g;
    s/\.rlen/\.query_alignment_length/g;
    s/([[\].0-9a-zA-Z]*)\.qqual/pysam.toQualityString($1\.query_alignment_qualities)/g;
    s/\.qstart/\.query_alignment_start/g;
    s/\.qend/\.query_alignment_end/g;
    s/\.qlen/\.query_alignment_length/g;
    s/\.mrnm/\.next_reference_id/g;
    s/\.rnext/\.next_reference_id/g;
    s/\.mpos/\.next_reference_start/g;
    s/\.rname/\.reference_id/g;
    s/\.isize/\.query_length/g;
    s/\.cigar/\.cigartuples/g unless (/\.cigarstring/);

    s/\.blocks/\.getBlocks()/g;
    s/\.aligned_pairs/\.getAlignedPairs()/g;
    s/\.inferred_length/\.getInferredQueryLength()/g;

    s/\.overlap()/\.getOverlap()/g;

    # PileupProxy
    s/\.n([^a-zA-Z])/\.nsegments$1/g;
    
    # if (/\.mrnm/ || /\.rnext/ || /\.mpos/ || /\.rname/)
    # { 
    #     warn "Deprecated tag $& at line $.\n";
    # }
    print;
}
