import sys
import os

IS_PYTHON3 = sys.version_info[0] >= 3

if IS_PYTHON3:
    from itertools import zip_longest
    from urllib.request import urlopen
else:
    from itertools import izip as zip_longest
    from urllib2 import urlopen


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


def checkURL(url):
    '''return True if URL is available.

    A URL might not be available if it is the wrong URL
    or there is no connection to the URL.
    '''
    try:
        urlopen(url, timeout=1)
        return True
    except:
        return False

