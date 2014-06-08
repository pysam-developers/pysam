import sys
import os

IS_PYTHON3 = sys.version_info[0] >= 3

if IS_PYTHON3:
    from itertools import zip_longest
else:
    from itertools import izip as zip_longest


def checkBinaryEqual(filename1, filename2):
    '''return true if the two files are binary equal.'''
    if os.path.getsize(filename1) != os.path.getsize(filename2):
        return False

    infile1 = open(filename1, "rb")
    infile2 = open(filename2, "rb")

    def chariter(infile):
        while 1:
            c = infile.read(1)
            if c == b"":
                break
            yield c

    found = False
    for c1, c2 in zip_longest(chariter(infile1), chariter(infile2)):
        if c1 != c2:
            break
    else:
        found = True

    infile1.close()
    infile2.close()
    return found

