import math

DEFAULT_BACKGROUND = {"A": 0.25, "G": 0.25, "T": 0.25, "C": 0.25}

# note the documentation comments are similar to Numpy-style docstrings, but I didn't have the time to learn the docstrings properly

def weigh_matrix(sample, background = DEFAULT_BACKGROUND):
    """Creates a weigh matrix that includes pseudo frequencies

    Parameters
    ----------
    sample: Sample
        The sample to use to construct the weigh matrix
    background: dict[str, float], optional
        The background frequencies of nucleotides. Defaults to 0.25 for all four.
    
    Returns
    -------
    list
        Returns a list of dict[str, float]s where result[i][n] is log-odds-score of nucleotide n appearing at position i
    """

    counts = sample.nucleotide_counts()
    sample_size = sample.sample_size()
    result = []
    for i in range(0, sample.seq_len()):
        result += [{}]
        for n in ['A', "G", "T", "C"]:
            result[i][n] = math.log2((counts[i][n] + background[n]) / (sample_size + 1) / background[n])
    return result

def relative_individual_information(target, sample, background = DEFAULT_BACKGROUND):
    """Calculates the relative individual information of the target sequence 
    
    Parameters
    ----------
    target: str
        The sequence to calculate the relative individual information for
    sample: Sample
        The sample used to calculate the weigh matrix
    background: dict[str, float], optional
        The background frequencies of nucleotides. Defaults to 0.25 for all four.

        
    Returns
    -------
    float
        The relative individual information of the target sequence if background frequencies were provided, or just
        the individual information if they were not
    """

    wm = weigh_matrix(sample, background)
    result = 0
    for i in range(0, len(target)):
        result += wm[i][target[i]]
    return result

def blank_nucleotide_table(len):
    """Returns a blank table for counting nucleotides

    Parameters
    ----------
    len: int
        The length of the sequence you are creating a table for

    Returns
    -------
    list
        A list of length len where every element is the dict {"A": 0, "G": 0, "T": 0, "C": 0}
    """

    nucleotide_table = []
    for _ in range(0, len):
        nucleotide_table += [{"A": 0, "G": 0, "T": 0, "C": 0}]
    return nucleotide_table

"""interface Sample

Both ListSample and HMMSample are Samples, so they have many of the same methods. This section
documents their shared methods, so the same documentation does not need to be copied twice.

Methods
-------
sample_size(): int
    Returns the sample size
seq_len(): int
    Returns the length of the sequences in the sample
nucleotide_counts(): list
    Returns a list of dict[str][float]s, similar to blank_nucleotide_table(len), but filled in with
    accurate counts that can be used to create a weigh table

"""


class ListSample:
    """A Sample constructed from a list of strings, where each string is a sequence. Implements the Sample interface

    Methods
    _______
    to_hmm(): dict[int][(list, int)]
        Returns an HMM encoding the sequences in this sample
    to_hmm_sample(): HMMSample
        Returns an HMMSample constructed from the hmm produced by calling to_hmm()
    (see Sample documentation for Sample methods)    
    """
    
    def __init__(self, sequences, start = 0, end = -1):
        """Initializes a ListSample, which implements the Sample interface

        This initializer takes optional start and end arguments. If they are included, the object will treat the 
        sequences as only existing in the range from start (inclusive) to end (exclusive). All nucleotides outside
        that range will be ignored by all methods of this class.

        Parameters
        ----------
        sequences: list
            A list of strings each representing DNA sequences. They should be aligned and all be the same length.
        start: int, optional
            The start position of the sequences to observe (inclusive). Defaults to 0.
        end: int, optional
            The end position of the sequences to observe (exclusive). If this parameter is -1, then it will
            be ignored and the sequences will be considered to the end. Defaults to -1.
        """

        self._sequences = sequences
        self._start = start
        if end == -1:
            self._end = len(sequences[0])
        else:
            self._end = end

    def sample_size(self):
        return len(self._sequences)
    
    def seq_len(self):
        return self._end - self._start

    def nucleotide_counts(self):        
        counts = blank_nucleotide_table(self.seq_len())
        for seq in self._sequences:
            for i in range(0, self.seq_len()):
                counts[i][seq[i + self._start]] += 1
        return counts
    
    def to_hmm(self):
        """Converts the sequence list into an HMM
        
        Only encodes the parts between the start and end

        Returns
        -------
        dict[int][(list, int)]
            An HMM
        """

        counts = self.nucleotide_counts()
        sample_size = self.sample_size()

        hmm = {}
        for i in range(0, self.seq_len()):
            entry = []
            for n in counts[i]:
                entry += [(n, counts[i][n] / sample_size)]
            next_index = -1
            if i != self.seq_len() - 1:
                next_index = i + 1
            hmm[i] = (entry, next_index)

        return hmm
    
    def to_hmm_sample(self):
        """Converts the ListSample to an HMMSample

        Specifically, returns HMMSample(self.to_hmm(), self.sample_size())

        Returns
        -------
        This sample converted to an HMMSample        
        """
        
        return HMMSample(self.to_hmm(), self.sample_size())
    
class HMMSample:
    def __init__(self, hmm, sample_size, start = 0, end = -1):
        """Initialized an HMMSample, which implements the Sample interface

        This initializer takes optional start and end arguments. If they are included, the object will treat the 
        sequences as only existing in the range from start (inclusive) to end (exclusive). All nucleotides outside
        that range will be ignored by all methods of this class.

        Parameters
        ----------
        hmm: dict[int][(list, int)]
            A list of strings each representing DNA sequences. They should be aligned and all be the same length.
        sample_size: int
            The size of the sample from which the HMM was generated. This is needed to accurately calculate
            pseudo-frequencies
        start: int, optional
            The start position of the sequences to observe (inclusive). Defaults to 0.
        end: int, optional
            The end position of the sequences to observe (exclusive). If this parameter is -1, then it will
            be ignored and the sequences will be considered to the end. Defaults to -1.
        """
        
        self._hmm = hmm
        self._sample_size = sample_size
        self._start = start
        if end == -1:
            self._end = len(hmm)
        else:
            self._end = end

    def sample_size(self):
        return self._sample_size
    
    def seq_len(self):
        return self._end - self._start
    
    def nucleotide_counts(self):
        counts = blank_nucleotide_table(self.seq_len())
        for i in range(0, self.seq_len()):
            for nucleotide, probability in self._hmm[i + self._start][0]:
                counts[i][nucleotide] = probability * self._sample_size
        return counts    


test_sequences = ['ACGTACGA', 'ACGTACGT', 'ACGTCCGT', 'ACGTCCAT', 
                  'ACAGGCAT', 'ACAGGCTT', 'ACAGTCTT', 'ACAGTCTT']
list_sample = ListSample(test_sequences)
print(relative_individual_information(test_sequences[0], list_sample))

print(relative_individual_information('AAAAAAAA', list_sample))

test_hmm = list_sample.to_hmm()

for i in test_hmm:
    print(f'hmm[{i}] = {test_hmm[i]}')

print(relative_individual_information(test_sequences[0], HMMSample(test_hmm, list_sample.sample_size())))