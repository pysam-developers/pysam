import os
import sys
import subprocess
import threading
import errno
import unittest

from pysam import AlignmentFile

IS_PYTHON2 = sys.version_info[0] == 2

DATADIR = os.path.abspath(os.path.join(
    os.path.dirname(__file__),
    "pysam_data"))


def alignmentfile_writer_thread(infile, outfile):
    def _writer_thread(infile, outfile):
        """read from infile and write to outfile"""
        try:
            i = 0
            for record in infile:
                outfile.write(record)
                i += 1
        except IOError as e:
            if e.errno != errno.EPIPE:
                pass
        finally:
            outfile.close()

    writer = threading.Thread(target=_writer_thread, args=(infile, outfile))
    writer.daemon = True
    writer.start()
    return writer


class StreamTest(unittest.TestCase):

    def stream_process(self, proc, in_stream, out_stream, writer):

        with AlignmentFile(proc.stdout) as infile:
            read = 0
            for record in infile:
                read += 1
        return 0, read

    @unittest.skipIf(IS_PYTHON2, "no context manager in py2")
    def test_text_processing(self):

        with subprocess.Popen('head -n200',
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              shell=True) as proc:

            in_stream = AlignmentFile('pysam_data/ex1.bam')
            out_stream = AlignmentFile(proc.stdin, 'wh', header=in_stream.header)
            writer = alignmentfile_writer_thread(in_stream,
                                                 out_stream)

            written, read = self.stream_process(proc,
                                                in_stream,
                                                out_stream,
                                                writer)
            self.assertEqual(read, 198)

    @unittest.skip("test contains bug")
    def test_samtools_processing(self):
        
        # The following test causes the suite to hang
        # as the stream_processor raises:
        # ValueError: file has no sequences defined (mode='r') - is it SAM/BAM format?
        # The whole setup then hangs during exception handling.
        with subprocess.Popen('samtools view -b -f 4',
                              stdin=subprocess.PIPE,
                              stdout=subprocess.PIPE,
                              shell=True) as proc:
        
            in_stream = AlignmentFile('pysam_data/ex1.bam')
            out_stream = AlignmentFile(proc.stdin, 'wb', header=in_stream.header)
            writer = alignmentfile_writer_thread(in_stream,
                                                 out_stream)

            written, read = self.stream_process(proc,
                                                in_stream,
                                                out_stream,
                                                writer)
            self.assertEqual(read, 35)
        

if __name__ == "__main__":
    unittest.main()
