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
DEFAULT_BACKGROUND = {"A": 0.25, "G": 0.25, "T": 0.25, "C": 0.25}

test_sequences = ['ACGTACGA', 'ACGTACGT', 'ACGTCCGT', 'ACGTCCAT', 
                  'ACAGGCAT', 'ACAGGCTT', 'ACAGTCTT', 'ACAGTCTT']

def weigh_matrix(sequences, start, end, background = DEFAULT_BACKGROUND):
    nucleotide_table = []
    for _ in range(start, end):
        nucleotide_table += [{"A": 0, "G": 0, "T": 0, "C": 0}]
    for seq in sequences:
        for i in range(start, end):
            nucleotide_table[i - start][seq[i]] += 1
    result = []
    for i in range(0, end - start):
        result += [{}]
        for n in ['A', "G", "T", "C"]:
            result[i][n] = math.log2((nucleotide_table[i][n] + background[n]) / (len(sequences) + 1) / background[n])
    return result


def relative_individual_information(target, sequences, start, end, background = DEFAULT_BACKGROUND):
    wm = weigh_matrix(sequences, start, end, background)
    result = 0
    for i in range(0, end-start):
        result += wm[i][target[i + start]]
    return result

print(relative_individual_information('ACGTACGA', test_sequences, 0, 8))