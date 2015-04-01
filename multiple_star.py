import numpy

# Global alignment dengan Needleman-Wunsch
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

# Menjajarkan seluruh pasangan pairwise
def all_pairs(s, sim_matrix, gap_penalty):

    score_matrix = numpy.zeros((len(s), len(s)), int)
    alignments = dict()

    for i in range(0, len(s) - 1):
        for j in range(i + 1, len(s)):
            score, alseq1, alseq2, _ = nw_align(s[i], s[j], sim_matrix, gap_penalty)
            score_matrix[i][j] = score
            score_matrix[j][i] = score

            alignments[(i, j)] = [alseq1, alseq2]
            alignments[(j, i)] = [alseq2, alseq1]

    # Pencarian center star
    sum_rows = score_matrix.sum(0)
    star_index = sum_rows.argmax()

    return star_index, alignments

# Fungsi tambahan untuk insert gap pada posisi tertentu
def insert_gap(s, index):
    return s[:index] + '-' + s[index:]

# Menggabungkan beberapa pairwise alignment menjadi satu multiple alignment
def merge_alignments(star_index, aligns):
    merged = []
    for k in aligns:
        if k[0] == star_index:
            if merged == []:
                merged.append(aligns[k][0])
                merged.append(aligns[k][1])
            else:
                star = merged[0]
                curr = aligns[k][0]
                for i in range(max(len(curr), len(star))):
                    if curr[i] == '-': # gap di curr, maka seluruh merged disisipkan gap
                        for m in range(len(merged)):
                            merged[m] = insert_gap(merged[m], i)
                        star = merged[0]
                    elif star[i] == '-': # gap di star, maka seluruh curr disisipkan gap
                        aligns[k][0] = insert_gap(aligns[k][0], i)
                        aligns[k][1] = insert_gap(aligns[k][1], i)
                        curr = aligns[k][0]
                    elif curr[i] == star[i]:
                        continue
                merged.append(aligns[k][1])

    return merged

if __name__ == "__main__":

    seq = ["ATTGCCATT", "ATGGCCATT", "ATCCAATTTT", "ATCTTCTT", "ATTGCCGATT"]

    sim_matrix = {
        ('A','A'): +2,
        ('G','A'): -1,  ('G','G'): +2,
        ('C','A'): -1,  ('C','G'): -1,  ('C','C'): +2,
        ('T','A'): -1,  ('T','G'): -1,  ('T','C'): -1,  ('T','T'): +2
    }
    gap_penalty = -1

    star, aligns = all_pairs(seq, sim_matrix, gap_penalty)

    merged = merge_alignments(star, aligns)

    for m in merged:
        print m











