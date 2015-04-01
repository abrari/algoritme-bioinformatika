import numpy

# Penjajaran global dengan algoritme Needleman-Wunsch
# Linear gap scoring, single path traceback
# sim_matrix adalah dictionary yang berisi subsitution matriks
def nw_align(seq1, seq2, sim_matrix, gap_penalty):
    rows, cols = len(seq1)+1, len(seq2)+1
    M = numpy.zeros((rows, cols), int)  # Matriks skor
    P = numpy.zeros((rows, cols), int)  # Matriks pointer traceback

    # Dynamic programming
    for i in range(rows):
        M[i][0] = i * gap_penalty
        P[i][0] = 3
    for j in range(cols):
        M[0][j] = j * gap_penalty
        P[0][j] = 2
    for i in range(1, rows):
        for j in range(1, cols):
            k1, k2 = seq1[i-1], seq2[j-1]
            if (k1, k2) in sim_matrix:
                d = M[i-1][j-1] + sim_matrix[(k1, k2)]
            else:
                d = M[i-1][j-1] + sim_matrix[(k2, k1)]
            u = M[i-1][j] + gap_penalty
            l = M[i][j-1] + gap_penalty
            M[i][j] = max(d, u, l)

            if M[i][j] == d: P[i][j] = 1
            elif M[i][j] == l: P[i][j] = 2
            elif M[i][j] == u: P[i][j] = 3

    # Traceback
    i, j = rows-1, cols-1
    al_seq1 = al_seq2 = ""
    al_score = M[i][j]
    while i > 0 and j > 0:
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

    while i > 0:
        al_seq1 = seq1[i-1] + al_seq1
        al_seq2 = "-" + al_seq2
        i -= 1

    while j > 0:
        al_seq1 = "-" + al_seq1
        al_seq2 = seq2[j-1] + al_seq2
        j -= 1

    al_symbol = ''
    for i in range(len(al_seq1)):
        al_symbol += "|" if al_seq1[i] == al_seq2[i] else " "

    # print M
    # print P

    return al_score, al_seq1, al_seq2, al_symbol

if __name__ == '__main__':
    seq1 = "ATGCTAGCTAGCTAGCTAGCTAGCATC"
    seq2 = "ATGCGGCATCCAGGGACTACTGATC"

    sim_matrix = {
        ('A','A'): +1,
        ('G','A'): -1,  ('G','G'): +1,
        ('C','A'): -1,  ('C','G'): -1,  ('C','C'): +1,
        ('T','A'): -1,  ('T','G'): -1,  ('T','C'): -1,  ('T','T'): +1
    }

    gap_penalty = -1

    score, al_seq1, al_seq2, al_sym = nw_align(seq1, seq2, sim_matrix, gap_penalty)

    print al_seq1
    print al_sym
    print al_seq2
    print "Score =", score
