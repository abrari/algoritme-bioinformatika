import re
from Bio import SeqIO


def prosite_to_regex(pattern):
    regex = pattern\
        .replace('{', '[^')\
        .replace('}', ']')\
        .replace('(', '{')\
        .replace(')', '}')\
        .replace(' ', '')\
        .replace('-', '')\
        .replace('x', '.')\
        .replace('>', '$')\
        .replace('<', '^')
    return regex

PS00014 = re.compile(prosite_to_regex('[KRHQSA]-[DENQ]-E-L>'))
p_falc_fasta = file('data/P_falc_RefProts.fa')

for record in SeqIO.parse(p_falc_fasta, 'fasta'):
    seq = str(record.seq)
    match = PS00014.search(seq)
    if match:
        print record.description
        a = match.start()
        b = match.end()
        print '%s at (%d,%d)' % (seq[a:b], a, b)
        print ''

PS00348 = re.compile(prosite_to_regex('M-C-N-S-S-C-[MV]-G-G-M-N-R-R'))
p53_fasta = file('data/p53_tumor.fa')

for record in SeqIO.parse(p53_fasta, 'fasta'):
    seq = str(record.seq)
    match = PS00348.search(seq)
    if match:
        print record.description
        a = match.start()
        b = match.end()
        print '%s at (%d,%d)' % (seq[a:b], a, b)
        print ''