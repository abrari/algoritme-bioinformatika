# -*- coding: utf-8 -*-

import numpy

# Penjajaran lokal dengan algoritme Smith-Waterman
# Linear gap scoring, multi path traceback
# sim_matrix adalah dictionary yang berisi subsitution matriks
def sw_align(seq1, seq2, sim_matrix, gap_penalty):
    rows, cols = len(seq1)+1, len(seq2)+1
    M = numpy.zeros((rows, cols), int)
    P = numpy.zeros((rows, cols), int) # 1 = diag, 2 = left, 3 = up

    # Dynamic programming
    for i in range(1, rows):
        for j in range(1, cols):
            k1, k2 = seq1[i-1], seq2[j-1]
            if (k1, k2) in sim_matrix:
                d = M[i-1][j-1] + sim_matrix[(k1, k2)]
            else:
                d = M[i-1][j-1] + sim_matrix[(k2, k1)]
            u = M[i-1][j] + gap_penalty
            l = M[i][j-1] + gap_penalty
            M[i][j] = max(d, u, l, 0)

            if M[i][j] == d: P[i][j] = 1
            elif M[i][j] == l: P[i][j] = 2
            elif M[i][j] == u: P[i][j] = 3

    # Traceback
    max_idx = numpy.where(M == M.max())     # cari indeks maksimum
    idx = zip(max_idx[0], max_idx[1])

    result = []

    for h in idx:   # untuk setiap posisi mulai trackback
        i, j = h
        al_seq1 = al_seq2 = ''
        score = M[i][j]
        while P[i][j] != 0:
            if P[i][j] == 1: # arah diagonal
                al_seq1 = seq1[i-1] + al_seq1
                al_seq2 = seq2[j-1] + al_seq2
                i -= 1
                j -= 1
            elif P[i][j] == 2: # arah kiri
                al_seq1 = '-' + al_seq1
                al_seq2 = seq2[j-1] + al_seq2
                j -= 1
            elif P[i][j] == 3:
                al_seq1 = seq1[i-1] + al_seq1
                al_seq2 = "-" + al_seq2
                i -= 1
            else:
                break

        al_symbol = ''
        for i in range(len(al_seq1)):
            al_symbol += "|" if al_seq1[i] == al_seq2[i] else " "

        result.append([score, al_seq1, al_seq2, al_symbol])

    return result, M, P

# Cetak matriks DP
def print_dp_matrix(seq1, seq2, M, P):
    rows, cols = len(seq1)+1, len(seq2)+1
    arrow = [' ', '↖', '←', '↑']

    seq2 = '-' + seq2
    print ' ',
    for j in range(0, cols):
        print ' ' + seq2[j],
    print

    seq1 = '-' + seq1
    for i in range(0, rows):
        print seq1[i],
        for j in range(0, cols):
            print arrow[P[i][j]] + str(M[i][j]),
        print

if __name__ == '__main__':
    seq1 = "ACCTAAGG"
    seq2 = "GGCTCAATCA"

    sim_matrix = {
        ('A','A'): +2,
        ('G','A'): -1,  ('G','G'): +2,
        ('C','A'): -1,  ('C','G'): -1,  ('C','C'): +2,
        ('T','A'): -1,  ('T','G'): -1,  ('T','C'): -1,  ('T','T'): +2
    }

    gap_penalty = -3

    alns, M, P = sw_align(seq1, seq2, sim_matrix, gap_penalty)

    print_dp_matrix(seq1, seq2, M, P)
    print

    for aln in alns:
        print aln[1]
        print aln[3]
        print aln[2]
        print "Score =", aln[0]
        print

