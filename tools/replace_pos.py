#!/usr/bin/env python
# Usage:
# replace_pos.py mecab_command JUMAN|IPA [mecab_options] < in.KNP > out.JDP
# replace_pos.py juman_command                           < in.KNP > out.JDP
#   replace gold POS in Kyoto Univ. Text Corpus w/ auto POS given by MeCab/JUMAN
#   * you may want to set mecab_optiuons to '-d your_target_dict'
import sys, re, subprocess

def set_charset (data): # input coding to tagger
    for codec in ['shift_jis','utf-8','euc_jp','iso2022-jp']:
        try:
            data.decode (codec)
            return codec
        except:
            continue;
    else:
        sys.exit ("replace_pos.py: cannot decide input coding.")

class Binfo:
    """ bunsetsu infomation """
    def __init__ (self, offset, head_id, ptype):
        self.offset, self.head_id, self.ptype = offset, head_id, ptype;
        self.morphemes = []
    def header (self): return "%d%s" % (self.head_id, self.ptype)
    def morph  (self): return '\n'.join (morph for morph in self.morphemes)

if len (sys.argv) < 2:
    sys.exit ("""Usage:
  python %s mecab_command JUMAN|IPA [mecab_options] < in.KNP > out.JDP
  python %s juman_command                           < in.KNP > out.JDP"""
              % (sys.argv[0], sys.argv[0]))

tagger = len (sys.argv) == 2 and "JUMAN" or "MeCab"
posset = "JUMAN"
if tagger == "JUMAN":
    sys.argv.append ("-b")
else:
    posset = sys.argv.pop (2)
    if posset != "JUMAN" and posset != "IPA":
        sys.exit ("Unknown POS set: %s\n" % posset)

tagger_charset = set_charset (subprocess.Popen ('echo "X" | %s' % ' '.join (sys.argv[1:]), shell=True, stdout=subprocess.PIPE).communicate ()[0])
t = subprocess.Popen (sys.argv[1:], stdin=subprocess.PIPE, stdout=subprocess.PIPE)
sent   = ''
text   = ''
binfo  = []
stat   = { 'success': 0, 'failed': 0 }
for line_ in sys.stdin:
    text += line_
    if line_[0] == 'E': # EOS
        try:
            if tagger_charset != 'euc_jp':
                text = text.decode ('euc_jp').encode (tagger_charset)
        except:
            stat['failed'] += 1
            sys.stderr.write ("failed to decode: %s\n" % text.split ('\n')[0])
            text = ''
            continue
        lines = text[:-1].split ('\n')
        header, footer = lines[0], lines[-1]
        for line in lines[1:-1]:
            if line[0] == '*':
                head_type = line.split (' ')[2]
                binfo.append (Binfo (len (sent), int (head_type[:-1]), head_type[-1]))
            else:
                morph = line.split(' ')[0]
                sent += morph
        binfo.append (Binfo (len (sent), -1, 'D')) # dummy
        offset = 0
        i = 0
        t.stdin.write (sent + "\n")
        morphs = []
        while True:
            morphs.append (t.stdout.readline ()[:-1])
            if morphs[-1] == "EOS":
                break
        for morph in morphs:
            surface = ""
            feature = ""
            m = re.search ("^(\S+)[\s|\t](.+)$", morph)
            if m:
                surface = m.group (1)
                feature = m.group (2)
            # unknown words lacks fields related to surface form
            # modify here to be consistent with the MeCab output
            if binfo[i].offset == offset: # correct bunsetsu chunk
                i += 1
            elif i + 1 < len (binfo):
                if binfo[i + 1].offset == offset and binfo[i - 1].head_id == i:
                    # assuming compounds; you may want to add more conditions
                    binfo[i-1].head_id = binfo[i].head_id
                    del (binfo[i])
                    for bi in filter (lambda bj: bj.head_id >= i, binfo):
                        bi.head_id -= 1
                    i += 1
                elif binfo[i + 1].offset < offset: # mismatched
                    # sys.stderr.write ("failed to convert: %s" % header)
                    stat['failed'] += 1
                    break
            if surface:
                if posset == "IPA" and len (feature.split (',')) != 9:
                    morph += ",*,*" # complement empty info for unknown words
                binfo[i-1].morphemes.append (morph)
                offset += len (surface)
        else:
            text = header + "\n"
            for j in range (len (binfo) - 1):
                text += "* %d %s\n" % (j, binfo[j].header ())
                text += "%s\n" % binfo[j].morph ()
            text += footer + "\n"
            try:
                if tagger_charset != 'euc_jp':
                    text = text.decode (tagger_charset).encode ('euc_jp')
            except:
                stat['failed'] += 1
                sys.stderr.write ("failed to decode: %s\n" % header)
            else:
                stat['success'] += 1
                sys.stdout.write (text)
        sent = ''
        text = ''
        binfo [:] = []

sys.stderr.write ("(success, failed) = (%d, %d)\n"
                  % (stat['success'], stat['failed']))
