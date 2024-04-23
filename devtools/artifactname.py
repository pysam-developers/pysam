#!/usr/bin/env python3

import os
import re
import sys

pattern = re.compile(r'-cp([^-]*)-cp[^-]*-([^_]*)[0-9_]*_([^.]*)')
vers = set()
plats = set()
arches = set()

for fname in sys.argv[1:]:
    m = pattern.search(fname)
    vers.add(int(m[1]))
    plats.add(m[2])
    arches.add(m[3])

plat = '-'.join(sorted(plats))
arch = '-'.join(sorted(arches))
ver  = '-'.join(map(str, sorted(vers)))

tag = os.environ.get('NAMETAG', 'none')
tag = f'-{tag}' if tag != 'none' else ''

print(f'artifactname=wheels-{plat}-{arch}-{ver}{tag}')
