import math

# a function to calculate a weigh matrix for a list of sequences
# the start (inclusive) and end(exclusive) are 0-indexed positions from
#   in each sequence to include
# the sequences should be aligned already
# the output is a list of tuples. Each element of the outer list corresponds
#   to a sequence position, the inner tuple is the weigh values for A, G, T, and C
#   respectively.
# Pseudocounts are included
# If you pass in background nucleotide frequencies, it calculates the weigh matrix
#   to be used for relative individual information. Otherwise, it calculates the weigh
#   matrix for regular individual information
def weigh_matrix(sequences, start, end, background = {"A": 0.25, "G": 0.25, "T": 0.25, "C": 0.25}):
    nucleotide_table = []
    for _ in range(start, end):
        nucleotide_table += [{"A": 1, "G": 1, "T": 1, "C": 1}]
    for seq in sequences:
        for i in range(start, end):
            nucleotide_table[i - start][seq[i]] += 1
    result = []
    for i in range(0, end - start):
        result[i] = {}
        for n in ['A', "G", "T", "C"]:
            result[i][n] = math.log2(nucleotide_table[i][n] / len(sequences) / background[n])
    return result


def relative_individual_information(sequences, start, end, background = {"A": 0.25, "G": 0.25, "T": 0.25, "C": 0.25}):
    wm = weigh_matrix(sequences, start, end, background)
    result = []
    for seq in sequences:
        rii = 0
        for i in range(0, end - start):
            rii += wm[i][seq[i + start]]
        result += [rii]
    return result