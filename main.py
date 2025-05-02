
# a function to calculate a weigh matrix for a list of sequences
# the start (inclusive) and end(exclusive) are 0-indexed positions from
#   in each sequence to include
# the sequences should be aligned already
# the output is a list of tupples. Each element of the outer list coresponds
#   to a sequence position, the inner tuple is the weigh values for A, G, T, and C
#   respectively.
# Pseudocounts are included
def weigh_matrix(sequences, start, end):
    nucleotide_table = []
    for _ in range(start, end):
        nucleotide_table += [{"A": 1, "G": 1, "T": 1, "C": 1}]
    for seq in sequences:
        for i in range(start, end):
            nucleotide_table[i - start][seq[i]] += 1
    