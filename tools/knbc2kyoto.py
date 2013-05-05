#!/usr/bin/env/python
# Usage: knbc2kyoto.py [format]
#   convert KNP output to JDP for KyotoCorpus (if arg2 given) format
import sys, re

header, data, head, body, t2c = "", "", [], [], {}
for line in sys.stdin:
    if line[0] == '#':
        header = line
    elif line[0] == 'E':
        try:
            head = [re.sub (r'^\d+', lambda m: t2c[m.group ()], h) for h in head]
            for i in xrange (len (body)):
                data += "* %d %s\n%s" % (i, head[i], body[i])
            sys.stdout.write ("%s%s%s" % (header, data, line))
        except KeyError:
            sys.stderr.write ("fail to convert: %s" % header)
        header, data, head, body, t2c = "", "", [], [], {}
    elif line[0] == '*' or line[0] == '+':
        if line[0] == '*':
            body[len (body):], head[len (head):] = [""], [""]
        head[-1] = line[:-1].split(' ')[1]
        t2c[str (len (t2c))] = str (len (head) - 1)
    else:
        f = re.split (r'\s+', line)[:11]
        if len (sys.argv) >= 2: # KNP
            if f[0] == f[2]:
                f[2] = '*'
            body[-1] += ' '.join (f[0:3] + f[3:10:2] + ['*']) + "\n"
        else:
            body[-1] += "%s\t%s\n" % (f[0], ','.join (f[3:10:2] + [f[2], f[1], '*']))
