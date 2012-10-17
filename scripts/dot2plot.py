#!/usr/bin/env python

import sys
import numpy as np
import matplotlib.pyplot as plt
from operator import itemgetter


def parse_dot_file(filename):
    results, confluences = extract_nodes(filename)
    connections = extract_connection(filename, results)
    return results, connections, confluences

def extract_nodes(filename):
    confluences = {}
    line_arrays = {}
    for i, line in enumerate(open(filename)):
        if not 'label' in line:
            continue

        ary = line.split('\"')[1].split('\\n')

        if len(ary) < 3:
            continue
        if  not 'phonon_relax' in ary[1]:
            continue

        tid = int(ary[0].split()[0].strip('[]'))

        if len(ary[0].split()) > 3:
            # "[348] confluence with [80]"
            tid_parent = int(ary[0].split()[3].strip('[]'))
            confluences[tid] = tid_parent
        else:
            tid_parent = None

        # ['[687] done', 'phonon_relax-9', 'P2/c --> C2/m', '-9.390433/6']                
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

    return sorted(results, key=itemgetter(0)), confluences

def extract_connection(filename, results):
    tids = [x[0] for x in results]
    connections = []
    for i, line in enumerate(open(filename)):
        ary = line.split()
        if len(ary) != 4:
            continue
        if ary[1] == '->':
            left = int(ary[0][1:])
            right = int(ary[2][1:])
            if left in tids:
                connections.append([left, right])
    return connections

def draw_point(x, y, tid, t, confluence=None):
    plt.plot(x, y, 'ro')
    if confluence:
        plt.text(x, y, t + " [%d=%d]" % (tid, confluence),
                 rotation=60, ha='left', va='bottom')
    else:
        plt.text(x, y, t + " [%d]" % tid, rotation=60, ha='left', va='bottom')

def draw_line(x, y, style):
    plt.plot(x, y, style)

def draw_line2(x, y, style, t):
    plt.plot([x[0], x[1]], [y[0], y[0]], style)
    plt.plot([x[1], x[1]], [y[0], y[1]], style)
    plt.text(x[1], (y[0] + y[1]) / 2, t, color='0.25', size='small', ha='center')

def collect_points_lines(points, lines, tid_pos, tid, pos,
                         energies, texts, tids,
                         connections, ig_tids):

    if not tid in tids or tid in ig_tids:
        return pos
    
    pos += 1
    tid_pos[tid] = pos

    rights = []
    for (left, right) in connections:
        if left == tid:
            rights.append(right)
    
    i = tids.index(tid)
    points.append([pos, energies[i], tid, texts[i].split()[2]])

    for right in rights:
        if not right in tids:
            continue
        
        left = tid
        i_1 = tids.index(left)
        i_2 = tids.index(right)
        t_1 = texts[i_1].split()[2]
        t_2 = texts[i_2].split()[2]
        y_1 = energies[i_1]
        y_2 = energies[i_2]

        if abs(y_1 - y_2) < 1e-10 and t_1 == t_2:
            continue

        pos = collect_points_lines(points, lines, tid_pos, right, pos,
                                   energies, texts, tids,
                                   connections, ig_tids)

        if right in tid_pos:
            x_1 = tid_pos[left]
            x_2 = tid_pos[right]

            if abs(y_1 - y_2) < 1e-10 and t_1 == t_2:
                lines.append([[x_1, x_2], [y_1, y_2], 'r:', texts[i_2].split()[0]])
            else:
                lines.append([[x_1, x_2], [y_1, y_2], 'c-', texts[i_2].split()[0]])

    return pos

def plot(results, connections, confluences, reference_energy, ig_tids):
    energies = []
    texts = []
    tids = []
    
    for x in results:
        if None in x:
            continue
        tids.append(x[0])
        energies.append(x[1] / x[2] - reference_energy)
        texts.append(x[3])

    tid_pos = {}
    points = []
    lines = []
    pos = collect_points_lines(points, lines, tid_pos, tids[0], 0,
                               energies, texts, tids,
                               connections, ig_tids)

    for x, y, style, t in lines:
        draw_line2(x, y, style, t)
    
    for x, y, tid, text in points:
        if tid in confluences:
            draw_point(x, y, tid, text, confluences[tid])
        else:
            draw_point(x, y, tid, text)

    plt.xlim(0, pos + 1)
    emax = max(energies)
    emin = min(energies)
    de = emax - emin
    plt.ylim(emin - de * 0.1, emax + de * 0.1)
    


from optparse import OptionParser
parser = OptionParser()
parser.set_defaults(reference_energy=0.0,
                    ignore_tids=None,
                    output_filename=None,
                    ymax=None,
                    ymin=None)
parser.add_option("--reference", dest="reference_energy",
                  type="float",
                  help="Reference energy to be set as zero")
parser.add_option("-o", "--output", dest="output_filename",
                  action="store", type="string",
                  help="Output filename of PDF plot")
parser.add_option("--igtids", dest="ignore_tids", action="store",
                  type="string", help="Draw without the TIDs")
parser.add_option("--ymax", dest="ymax",
                  type="float",
                  help="Maximum y value of plot area")
parser.add_option("--ymin", dest="ymin",
                  type="float",
                  help="Minimum y value of plot area")
(options, args) = parser.parse_args()

if options.ignore_tids:
    ig_tids = [int(x) for x in options.ignore_tids.split()]
else:
    ig_tids = []

results, connections, confluences = parse_dot_file(args[0])
plot(results, connections, confluences, options.reference_energy, ig_tids)
plt.xticks([])
if options.ymax:
    plt.ylim(ymax=options.ymax)
if options.ymin:
    plt.ylim(ymin=options.ymin)
plt.grid(True)

if options.output_filename:
    plt.rcParams['backend'] = 'PDF'
    plt.rcParams['pdf.fonttype'] = 42
    plt.savefig(options.output_filename)
else:
    plt.show()
