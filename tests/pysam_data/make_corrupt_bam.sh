#!/bin/bash


# Generates a BAM with a specific type of corruption, where a BAM alignment record spans the
# boundary into a corrupt BGZF block. The other BAMs created in the makefile are generated in such a
# way that every alignment record is wholly contained within a single BGZF block; corruption of a
# BGZF block in this case is interpretted by htslib as a "truncation" (return code = -2), whereas
# the boundary into a corrupt BGZF block being in the middle of an alignment record causes htslib to
# return different error codes (< -2).


INPUT_BAM=$1
OUTPUT_BAM=$2

# find the start of the second BGZF block since the first block often has only header information;
# BGZF blocks in BAM files all start with this sequence of bytes: 0x1f 0x8b 0x08 0x04
OFFSET=$(grep -abo $'\x1f\x8b\x08\x04' $INPUT_BAM | grep -aoP '^\d+' | sed -n 2p)

# decompress the blocks starting at the offset and split off the first ~35kb
# into one file and the rest into a second (typical decompressed block size is
# about 64kb)
tail -c +$(($OFFSET + 1)) $INPUT_BAM | bgzip -dc > _chunk_decompressed
head -c 35000 < _chunk_decompressed > _chunk1
tail -c +35001 < _chunk_decompressed > _chunk2

# recompress and corrupt the second chunk with stray null bytes ~5k into the compressed file
bgzip -c < _chunk2 > _chunk2.bgzip
dd if=/dev/zero of=_chunk2.bgzip bs=1 seek=5000 count=10 conv=notrunc

# glue everything back together:
#  * the skipped header
#  * the uncorrupt 35kb chunk
#  * the remainder with the corrupt BGZF block
head -c $OFFSET $INPUT_BAM > $OUTPUT_BAM
bgzip -c < _chunk1 >> $OUTPUT_BAM
cat _chunk2.bgzip >> $OUTPUT_BAM

# clean up
rm _chunk*
