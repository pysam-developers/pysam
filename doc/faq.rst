===
FAQ
===

pysam coordinates are wrong
===========================

pysam uses 0-based coordinates and the half-open notation for ranges
as does python. Coordinates and intervals reported from pysam always
follow that convention.

Confusion might arise as different file formats might have different
conventions. For example, the SAM format is 1-based while the BAM
format is 0-based. It is important to remember that pysam will always
conform to the python convention and translate to/from the file format
automatically.

The only exception is the :term:`region` string in the
:meth:`~pysam.AlignmentFile.fetch` and
:meth:`~pysam.AlignmentFile.pileup` methods. This string follows the
convention of the samtools command line utilities. The same is true
for any coordinates passed to the samtools command utilities directly,
such as :meth:`pysam.mpileup`.

BAM files with a large number of reference sequences is slow
============================================================

If you have many reference sequences in a bam file, the following
might be slow::

      track = pysam.AlignmentFile(fname, "rb")
      for aln in track.fetch():
      	  pass
	  
The reason is that track.fetch() will iterate through the bam file
for each reference sequence in the order as it is defined in the
header. This might require a lot of jumping around in the file. To
avoid this, use::

      track = pysam.AlignmentFile(fname, "rb")
      for aln in track.fetch( until_eof = True ):
      	  pass
 
This will iterate through reads as they appear in the file.

Weirdness with spliced reads in samfile.pileup(chr,start,end) given spliced alignments from an RNA-seq bam file
===============================================================================================================

Spliced reads are reported within samfile.pileup. To ignore these
in your analysis, test the flags ``is_del == True and indel=0``
in the :class:`~.PileupRead` object.

I can't edit quality scores in place
====================================

Editing reads in-place generally works, though there is some
quirk to be aware of. Assigning to AlignedRead.seq will invalidate 
any quality scores in AlignedRead.qual. The reason is that samtools
manages the memory of the sequence and quality scores together 
and thus requires them to always be of the same length or 0.

Thus, to in-place edit the sequence and quality scores, copies of
the quality scores need to be taken. Consider trimming for example::

    q = read.qual
    read.seq = read.seq[5:10]
    read.qual = q[5:10]
 

Why is there no SNPCaller class anymore?
=========================================

SNP calling is highly complex and heavily parameterized. There was a
danger that the pysam implementations might show different behaviour from the
samtools implementation, which would have caused a lot of confusion.

The best way to use samtools SNP calling from python is to use the 
:meth:`pysam.mpileup` command and parse the output  directly.

I get an error 'PileupProxy accessed after iterator finished'
=============================================================

Pysam works by providing proxy objects to objects defined within
the C-samtools package. Thus, some attention must be paid at the
lifetime of objects. The following to code snippets will cause an
error::

    s = AlignmentFile('ex1.bam')
    for p in s.pileup('chr1', 1000,1010):
        pass
    
    for pp in p.pileups:
        print pp

The iteration has finished, thus the contents of p are invalid. A
variation of this::

    p = next(AlignmentFile('ex1.bam').pileup('chr1', 1000, 1010))
    for pp in p.pileups:
        print pp

Again, the iteration finishes as the temporary iterator created
by pileup goes out of scope. The solution is to keep a handle
to the iterator that remains alive::

    i = AlignmentFile('ex1.bam').pileup( 'chr1', 1000, 1010)
    p = next(i)
    for pp in p.pileups:
        print pp

Psyam won't compile
===================

Compiling pysam can be tricky as there are numerous variables that
differ between build environments such as OS, version, python version,
and compiler. It is difficult to build software that build cleanly
on all systems and the process might fail. Please see the 
`pysam user group
<https://groups.google.com/forum/#!forum/pysam-user-group>`_
for common issues.

If there is a build issue, read the generated output carefully -
generally the cause of the problem is among the first errors to be
reported. For example, you will need to have the development version
of python installed that includes the header files such as
:file:`Python.h`. If that file is missing, the compiler will report
this at the very top of its error messages but will follow it 
with any unknown function or variable definition it encounters later
on.

A general advice is to always use the latest version on python_ and
cython_ when building pysam. There are some known incompatibilities:

* Python 3.4 requires cython 0.20.2 or later (see `here
  <https://github.com/pysam-developers/pysam/issues/37>`_)




