#!/usr/bin/env python

import sys

def reformat_dot_file(filename):
    for i, line in enumerate(open(filename)):
        if not 'label' in line:
            continue
        arr = line.split('\"')[1].split('\\n')
        if len(arr) > 2:
            print i, arr[1], arr[2]


from optparse import OptionParser
parser = OptionParser()
parser.set_defaults(split=False)
parser.add_option("--split", dest="split",
                  action="store_true",
                  help="Split and reformat .dot file")
(options, args) = parser.parse_args()

if options.split:
    reformat_dot_file(args[0])
