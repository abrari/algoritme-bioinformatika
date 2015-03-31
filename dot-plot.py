import numpy

# Menghitung jumlah karakter yang sama di str1 dan str2
def count_match(str1, str2):
    match = 0
    for i in range(len(str1)):
        if str1[i] == str2[i]:
            match += 1
    return match

# Pembuatan dotplot dengan window dan theshold tertentu
def compute_dotplot(str1, str2, win, threshold):
    dot_plot = numpy.zeros((len(str2), len(str1)), int)

    for i in range(len(str2) - win + 1):
        for j in range(len(str1) - win + 1):

            # Jika similarity dari window mencapai threshold, maka jadi satu dot
            if count_match(str2[i:i+win], str1[j:j+win]) >= threshold:
                dot_plot[i][j] = 1

    return dot_plot

# Menggambar dotplot
def print_dotplot(str1, str2, dot_plot):
    row = ""
    for ch in str1:
        row = row + "" + ch + ""
    print ("    " + row)
    print ("  +" + "-"*len(str1) + "-")
    row = ""
    for i in range(len(str2)):
        row = row + str2[i] + " | "
        for j in range(len(str1)):
            dot = "#" if dot_plot[i][j] == 1 else " "
            row = row + dot
        print (row)
        row = ""

if __name__ == "__main__":
    seq1 = "TAGCTAGCTATCCATCCCCATCCACATGTAC"
    seq2 = "TAGCTAGCTATCCATCCATGGAC"

    window = 5
    threshold = 4

    plot = compute_dotplot(seq1, seq2, window, threshold)
    print_dotplot(seq1, seq2, plot)