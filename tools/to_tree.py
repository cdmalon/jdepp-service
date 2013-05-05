#!/usr/bin/env python
# -*- coding: utf-8 -*-
# Usage: to_tree.py [-p, -prob] < in.JDP
#   translate parser output into human-readable dependency tree structure
import re, sys, os

if "-h" in sys.argv[1:] or "--help" in sys.argv[1:]:
    sys.exit ("Usage: %s\n  -p, --prob\tshow dependency probability" % sys.argv[0])
prob = "-p" in sys.argv[1:] or "--prob" in sys.argv[1:]

# customizable parameters
indent = 4          # for one dependency arc
offset = 2          # offset from the left
color  = "\033[31m" # color (red) for incorrect dependency arc

def set_charset (data): # input coding
    for codec in ['shift_jis','utf-8','euc_jp','iso2022-jp']:
        try:
            data.decode (codec)
            return codec
        except:
            continue;
    else:
        sys.exit ("to_tree.py: cannot decide input coding.")

class Binfo:
    """ bunsetsu infomation """
    def __init__ (self, *args):
        self.id, self.head, self.prob, self.fail, self.gold = args
        self.morph, self.depth, self.first_child = "", 0, -1
    def len    (self)         : return len (unicode (self.morph, 'utf-8'))
    def offset (self, offset) : return offset - self.len () - self.depth

def treeify (binfo):
    tree = ""
    for c in reversed (binfo[:-1]):
        c.depth = binfo[c.head].depth + indent
        binfo[c.head].first_child = c.id
    width = offset + max (b.len () + b.depth for b in binfo) # tree width
    for b in binfo:
        if b.head == -1:
            if b.id != len (binfo) - 1:
                sys.exit ("no head information; chunking output [-I 1]?")
            tree += "% 3d:%s%s" % (b.id, "　" * b.offset (width), b.morph)
        else:
            tree += b.fail and color or ""
            tree += "% 3d:%s%s" % (b.id, "　" * b.offset (width), b.morph)
            h = binfo[b.head]
            tree += "━" * (b.depth - h.depth - (prob and b.prob > 0 and 3 or 1))
            tree += prob and b.prob > 0 and "%.2f" % b.prob or "" # "─"
            tree += b.id == h.first_child and "┓" or "┫" # "┐" or "┤"
            tree += b.fail and "%-4s\033[0m" % b.gold or ""
            while h.head != -1: # draw arcs spanning from x < b to y > h
                c, h = h, binfo[h.head]
                tree += "　" * (c.depth - h.depth - (b.fail and 3 or 1))
                tree += h.first_child < b.id and "┃" or "　" # "│" or "　" 
                b.fail = False
            tree += "\n"
    return tree

binfo   = []
text    = ""
charset = ''
for line in iter (sys.stdin.readline, ""): # no buffering
    text += line
    if line[0] == 'E': # EOS
        if not charset: # set charset
            charset = set_charset (text)
        if charset != 'utf-8':
            text = text.decode (charset).encode ('utf-8')
        lines = text[:-1].split ('\n')
        header, footer = lines[0], lines[-1]
        for line_ in lines[1:-1]:
            if line_[0] == '*':
                h, p, g = re.split (r'[\s|@]', line_ + "@-1 -1D")[2:5]
                binfo.append (Binfo (len (binfo), int (h[:-1]), float (p),
                                     g[:-1] != "-1" and h[:-1] != g[:-1], g))
            else:
                surface = re.split ('[\t ]', line_)[0]
                binfo[-1].morph += surface
        text = '\n'.join ([header, treeify (binfo) + footer])
        if charset != 'utf-8':
            text = text.decode ('utf-8').encode (charset)
        sys.stdout.write (text + "\n")
        binfo[:] = []
        text  = ""
