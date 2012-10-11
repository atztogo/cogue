#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt

def split_dot_file(filename):
    line_arrays = {}
    for i, line in enumerate(open(filename)):
        if not 'label' in line:
            continue
        ary = line.split('\"')[1].split('\\n')
        if len(ary) > 2:
            tid = int(ary[0].split()[0].strip('[]'))
            if len(ary[0].split()) > 3:
                tid_parent = int(ary[0].split()[3].strip('[]'))
            else:
                tid_parent = None
            if len(ary) == 4:
                energy = float(ary[3].split('/')[0])
                num_atom = int(ary[3].split('/')[1])
            else:
                energy = None
                num_atom = None

            line_arrays[tid] = [tid_parent, energy, num_atom, ary[2]]


    for tid in line_arrays:
        if line_arrays[tid][0]:
            line_arrays[tid][1] = line_arrays[line_arrays[tid][0]][1]
            line_arrays[tid][2] = line_arrays[line_arrays[tid][0]][2]

    results = []
    for k, v in line_arrays.iteritems():
        results.append([k, v[1], v[2], v[3]])

    return results

def plot(results):
    tids = [x[0] for x in results]
    energies = [(x[1] / x[2]) for x in results]
    texts = [x[3] for x in results]

    for i, (x, y, t) in enumerate(zip(tids, energies, texts)):
        plt.plot(i + 1, y, 'o')
        plt.text(i + 1, y, t.split()[2] + " [%d]" % x, rotation=45, ha='left', va='bottom')

    plt.xlim(0, len(tids) + 1)
    plt.xticks([])
    plt.grid(True)
    

from optparse import OptionParser
parser = OptionParser()
parser.set_defaults(split=False)
parser.add_option("--split", dest="split",
                  action="store_true",
                  help="Split and reformat .dot file")
(options, args) = parser.parse_args()

results = split_dot_file(args[0])
plot(results)
plt.show()
