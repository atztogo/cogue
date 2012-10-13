#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter

def extract_nodes(filename):
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
        # tid, energy, num_atom, space_group_string
        results.append([k, v[1], v[2], v[3]])

    return sorted(results, key=itemgetter(0))

def extract_connection(filename, results):
    connections = []
    for i, line in enumerate(open(filename)):
        ary = line.split()
        if len(ary) != 4:
            continue
        if ary[1] == '->':
            tids = [x[0] for x in results]
            left = int(ary[0][1:])
            right = int(ary[2][1:])
            if left in tids:
                connections.append([left, right])
    return connections
        
def parse_dot_file(filename):
    results = extract_nodes(filename)
    connections = extract_connection(filename, results)
    return results, connections

def plot(results, connections, reference_energy):
    tids = [x[0] for x in results]
    energies = [x[1] / x[2] - reference_energy for x in results]
    texts = [x[3] for x in results]
    tid_pos = {}
    points = []

    def draw_point(x, y, tid, t):
        plt.plot(x, y, 'ro')
        plt.text(x, y, t.split()[2] + " [%d]" % tid, rotation=45, ha='left', va='bottom')

    def recr(tid, pos):
        pos += 1
        tid_pos[tid] = pos

        rights = []
        for (left, right) in connections:
            if left == tid:
                rights.append(right)
        
        i = tids.index(tid)
        points.append([pos, energies[i], tid, texts[i]])

        for right in rights:
            pos = recr(right, pos)

        return pos

    recr(tids[0], 0)

    for (left, right) in connections:
        i_1 = tids.index(left)
        x_1 = tid_pos[left]
        y_1 = energies[i_1]
        i_2 = tids.index(right)
        x_2 = tid_pos[right]
        y_2 = energies[i_2]
        plt.plot([x_1, x_2], [y_1, y_2], 'c-')

    for x, y, tid, text in points:
        draw_point(x, y, tid, text)

    plt.xlim(0, len(results) + 1)
    emax = max(energies)
    emin = min(energies)
    de = emax - emin
    plt.ylim(emin - de * 0.1, emax + de * 0.1)
    

from optparse import OptionParser
parser = OptionParser()
parser.set_defaults(reference_energy=0.0,
                    output_filename=None)
parser.add_option("--reference", dest="reference_energy",
                  type="float",
                  help="Reference energy to be set as zero")
parser.add_option("-o", "--output", dest="output_filename",
                  action="store", type="string",
                  help="Output filename of PDF plot")
(options, args) = parser.parse_args()

results, connections = parse_dot_file(args[0])
plot(results, connections, options.reference_energy)
plt.xticks([])
plt.grid(True)

if options.output_filename:
    plt.rcParams['backend'] = 'PDF'
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig(options.output_filename)
else:
    plt.show()
