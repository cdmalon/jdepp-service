#!/usr/bin/env python -u
# -*- coding: utf-8 -*-
# Usage: to_chunk.py [-p, -prob] < in.JDP
#   translate chunker output into human-readable chunked sequences
import sys, re
if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    sys.exit ("Usage: %s\n  -p, --prob\tshow chunking probability" % sys.argv[0])
prob = "-p" in sys.argv[1:] or "--prob" in sys.argv[1:]

flag = False
bid  = 0
for line in iter (sys.stdin.readline, ""): # no buffering
    if line[0] == '#' or line[0] == 'E':
        sys.stdout.write (line)
        bid = 0
    elif line[0] == '*':
        flag = True
    else:
        field = re.split ('[\t ]', line[:-1])
        surf, (auto, gold) = field[0], field[-2:]
        prob &= auto[0] == 'B' or auto[0] == 'I'
        p = prob and "%.2f" % float (auto[2:]) or ""
        if bid > 0:
            if auto[0] == 'B' and gold == 'B':
                sys.stdout.write ("│%s " % p) # "┃"
            elif auto[0] == 'I' and gold == 'B': # false negative
                sys.stdout.write ("\033[34m│%s \033[0m" % p)
            elif auto[0] == 'B' and gold == 'I': # false positive
                sys.stdout.write ("\033[31m│%s \033[0m" % p)
            elif flag:
                sys.stdout.write ("│")
        bid += flag and 1 or 0
        flag = False
        sys.stdout.write ("%s " % surf)
