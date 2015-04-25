# Modified from https://github.com/martijnvermaat/pymotif

from gibbs import Gibbs
from Bio import SeqIO

if __name__ == '__main__':

    motif_length = 10
    dna_data = 'data/MA0005.1.fa'
    max_iter = 100

    sequences = [{'title':          record.description,
                  'sequence':       str(record.seq).upper(),
                  'motif_position': 0}
                 for record in SeqIO.parse(file(dna_data), 'fasta')]

    g = Gibbs(sequences, motif_length)

    g.find_motif(max_iter)

    print "Motif occurrences in sequences follow"

    for i in range(len(sequences)):
        start, end = (sequences[i]['motif_position'],
                      sequences[i]['motif_position'] + motif_length)
        print "Sequence #%2i  %s  (at position %i)" % (
            i + 1,
            sequences[i]['sequence'][start:end],
            sequences[i]['motif_position'] + 1)




