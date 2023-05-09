import numpy as np
import sys

# Define scoring system
def SubstitutionMatrix(x, y):
    if x == y:
        return 2
    else:
        return -1


gap_penalty = -1

# Needleman Wunsch algorithm
def NeedlemanWunsch(seq1, seq2):
    n, m = len(seq1), len(seq2)
    dpMatrix = np.zeros((n + 1, m + 1), dtype=int)

    for i in range(1, n + 1):
        dpMatrix[i, 0] = gap_penalty * i
    for j in range(1, m + 1):
        dpMatrix[0, j] = gap_penalty * j

    for i in range(1, n + 1):
        for j in range(1, m + 1):
            match = dpMatrix[i - 1, j - 1] + SubstitutionMatrix(seq1[i - 1], seq2[j - 1])
            delete = dpMatrix[i - 1, j] + gap_penalty
            insert = dpMatrix[i, j - 1] + gap_penalty
            dpMatrix[i, j] = max(match, delete, insert)

    alignedSeq1 = []
    alignedSeq2 = []
    i, j = n, m

    while i > 0 or j > 0:
        if i > 0 and j > 0 and dpMatrix[i, j] == dpMatrix[i - 1, j - 1] + SubstitutionMatrix(seq1[i - 1], seq2[j - 1]):
            alignedSeq1.append(seq1[i - 1])
            alignedSeq2.append(seq2[j - 1])
            i -= 1
            j -= 1
        elif i > 0 and dpMatrix[i, j] == dpMatrix[i - 1, j] + gap_penalty:
            alignedSeq1.append(seq1[i - 1])
            alignedSeq2.append('-')
            i -= 1
        else:
            alignedSeq1.append('-')
            alignedSeq2.append(seq2[j - 1])
            j -= 1

    alignedSeq1 = ''.join(alignedSeq1[::-1])
    alignedSeq2 = ''.join(alignedSeq2[::-1])

    return alignedSeq1, alignedSeq2, dpMatrix[n, m]

if __name__ == "__main__":
    seq1 = input("Enter sequence 1: ").upper()
    seq2 = input("Enter sequence 2: ").upper()

    aligned_seq1, aligned_seq2, score = NeedlemanWunsch(seq1, seq2)
    print("\nSequence 1:", aligned_seq1)
    print("Sequence 2:", aligned_seq2)
    print("Alignment score:", score)