# Needleman-Wunsch-Algorithm

Implementation of a DNA sequence alignment algorithm. Time and Space Complexity: O(nm) where n and m are the lengths of the sequences.

## Background:
Strands of DNA are made up of molecules called nucleotides, of which there are four types: Adenosine (A), Thymine (T), Cytosine (C), annd Guanine (G). Our genome is made up of a long sequence of millions of these nucleotides with a variable ordering. Some parts are ordered the same, making up specific genes that all of us possess. Other areas, such as those that code for aspects such as appearance, are variable even among the same species. The further two species are in similarity and genetic structure, the greater the difference in their DNA sequence. 
An important part of biological research/application is comparing two DNA sequences and assessing similarity. From an academic perspective, this can be used to determine the phylogeny of some unknown DNA sequence, or the relationship between two known species' DNA. It can also be used to show how a certain DNA sequence mutated (insertion, deletion, etc) to form a new one, and where the most likely changes occured along the chain. Assuming two sequences of size N and M (where size is measured in number of nucleotides), the time complexity of a comprehensive algorithm that checks every possible comparison of nucleotides at each index is O(NM) (or if the sequences are the same length, O(n^2). With three or more sequences, the problem is NP-Complete. If we assume the sequences are the same size N for simplicity, we would need a matrix of size N^N to check every combination of nucleotides at a given index with every other nucleotide of every other sequence. The time and space complexity is therefore O(N^N). The Needleman-Wunsch algorithm utilizes dynamic programming to compare two DNA sequences and line them up based on matching nucleotides.

## How It Works
The DNA strands of interest (for now, we'll work with just two) are lined up along the row and column of a matrix, leaving an empty space at the first cell of each, as follows:

(Using sample sequences ATTCCG and ATCG)

|   | o | A | T | T | C | C | G |
| - | - | - | - | - | - | - | - |
| o |   |   |   |   |   |   |   |
| A |   |   |   |   |   |   |   |
| T |   |   |   |   |   |   |   |
| C |   |   |   |   |   |   |   |
| G |   |   |   |   |   |   |   |

Values of the matrix are initialized with a 0 in the top leftmost corner, followed by -1, -2...-n along the first row and -1, -2....-m down the first column. These numbers represent the "score" at that comparison, which will be utilized later.


|   | o | A | T | T | C | C | G |
| - | - | - | - | - | - | - | - |
| o | 0 |-1 |-2 |-3 |-4 |-5 |-6 |
| A |-1 |   |   |   |   |   |   |
| T |-2 |   |   |   |   |   |   |
| C |-3 |   |   |   |   |   |   |
| G |-4 |   |   |   |   |   |   |

Next, at each cell, we compare the corresponding nucleotides at the row and column index. We'll use a scoring matrix assigned as follows.

Match = +1
Mismatch/Insertion/Deletion = -1

These values can be assigned arbitrarily (a match should always be better than a mismatch/indel). We can use a generic matrix like the one above or, in real-world situations, an expert biologist would assign one for us. Recall that we are comparing the similarity of two DNA sequences, either how likely they are to be genetically related or how viable it is that one is the ancestor of another. During reproduction, the child DNA mostly remains the same when copied from the parents, but in rare cases a mistake is made during this process. These rare cases are actually the basis for evolution and mutation (and why we aren't all clones of each other). If the nucleotide didn't change between sequences, we'll refer to it as a Match. In our example sequences, only the first two nucleotides match:

ATTCCG
ATCG

If the nucleotides don't match, then it was caused by either a mutation, and insertion, or a deletion. Either the wrong nucleotide was copied, an extra nucleotide was inserted, or a nucleotide that should've been copied was never added to the chain. Again, simple biological mistakes like this do occur, and they are very rare. If our hypothesis is that the first sequence is the ancestor, then we can see in this example that there were two deletions that occured. These will be represented as gaps in our final alignment. Mismatches will tell us where the most likely gaps occured, while maximizing the number of matches.

To show how to apply the scores, we'll start with the first empty cell in the top leftmost corner. Three scores are generated based on values in neighboring cells, and the max of the three is chosen:

|   | o | A | T | T | C | C | G |
| - | - | - | - | - | - | - | - |
| o | 0 |-1 |-2 |-3 |-4 |-5 |-6 |
| A |-1 | X |   |   |   |   |   |
| T |-2 |   |   |   |   |   |   |
| C |-3 |   |   |   |   |   |   |
| G |-4 |   |   |   |   |   |   |

Score 1: Top-left neighbor (0) + Match, A and A (1) = 1 <br />
Score 2: Top neighbor (-1) + Indel, A and empty cell (-1) = -2 <br />
Score 3: Left neighbor (-1) + Indel, A and empty cell (-1) = -2

Highest score is 1, so it's recorded in the matrix. An imaginary "arrow" is drawn to the cell where the score came from. This is used later during th traceback step.

|   | o | A | T | T | C | C | G |
| - | - | - | - | - | - | - | - |
| o | 0 |-1 |-2 |-3 |-4 |-5 |-6 |
| A |-1 | 1 | X |   |   |   |   |
| T |-2 |   |   |   |   |   |   |
| C |-3 |   |   |   |   |   |   |
| G |-4 |   |   |   |   |   |   |

We can apply the same method to the cell directly to the right:

Score 1: Top-left neighbor (-1) + Mismatch, T and A (-1) = -2 <br />
Score 2: Top neighbor (-2) + Indel (-1) = -3 <br />
Score 3: Left neighbor (1) + Indel (-1) = 0

Once again we update with the max of the three scores.

|   | o | A | T | T | C | C | G |
| - | - | - | - | - | - | - | - |
| o | 0 |-1 |-2 |-3 |-4 |-5 |-6 |
| A |-1 | 1 | 0 | X |   |   |   |
| T |-2 |   |   |   |   |   |   |
| C |-3 |   |   |   |   |   |   |
| G |-4 |   |   |   |   |   |   |

One more example, next cell on the right.

Score 1: Top-left neighbor (-2) + Mismatch, T and A (-1) = -3 <br />
Score 2: Top neighbor (-3) + Indel (-1) = -4 <br />
Score 3: Left neighbor (0) + Indel (-1) = -1

|   | o | A | T | T | C | C | G |
| - | - | - | - | - | - | - | - |
| o | 0 |-1 |-2 |-3 |-4 |-5 |-6 |
| A |-1 | 1 | 0 |-1 |   |   |   |
| T |-2 |   |   |   |   |   |   |
| C |-3 |   |   |   |   |   |   |
| G |-4 |   |   |   |   |   |   |

You may recognize a pattern, the next cell (for this row at least) depends on its left neighbor since the other
scores will always be lower. If we fill out the rest of the table, we get:


|   | o | A | T | T | C | C | G |
| - | - | - | - | - | - | - | - |
| o | 0 |-1 |-2 |-3 |-4 |-5 |-6 |
| A |-1 | 1 | 0 |-1 |-2 |-3 |-4 |
| T |-2 | 0 | 2 | 1 | 0 |-1 |-2 |
| C |-3 |-1 | 1 | 1 | 2 | 1 | 0 |
| G |-4 |-2 | 0 | 0 | 1 | 1 | 2 |

The next step is to traceback through the matrix, starting at the bottom rightmost cell. We  follow the arrows back until we reach the last cell. With a path set, start with the top leftmost cell, pointed to by the final arrow and form the alignment by travelling back down:

If we move diagonally, write both nucleotides. If we move down, there's a gap in the first (row) sequence. If we move right, there's a gap in the second (column) sequence. Applying this to the matrix above, we get the following alignment:

ATTCCG <br />
AT-C-G

As we see, all possible matching nucleotides were aligned, as optimial as possible with respect to the scoring system. With longer and more complex chains (with multiple matchings), general optimality will be determined by the scoring system (perhaps a match is given more weight than a mismatch or indel).

Given the very general nature of the algorithm (this code itself), it can be run on any seqeuence of text, including RNA, protien chains, or something else entirely. It'll be up to the user to infer the output based on their parameters, but the general principle remains the same. However, scoring systems accounting for DNA base pairs, which abide by basic chemical laws, will be designed much differently than amino acid chains or plain text.


## Running the Code

The code is written in Perl 5, one of the preferred languages for bioinformatics-related programming.

Download NWA.pl along with test1.txt and test2.txt (you can also create these files on your own).

```
c/user$ ls
NWA.pl
test1.txt
test2.txt
```

Populate test1 and test2 with the first and second DNA sequence (respectively) you want to compare. Include only the DNA sequence of nucleotides.

Here's a databse that allows you to search based on organism/gene: https://www.ncbi.nlm.nih.gov/genbank/ Search for a sequence, select the result, and click "FASTA" below the name. Copy this sequence (without label) into the test files.


In the cwd, run "perl NWA.pl test1.txt test2.txt". When the run is complete, the program will print the original sequence, the scoring matrix, the traceback matrix, and the alignment along with final score (positive usually means good).
 (In the traceback matrix, each letter in the cell tells you where the arrow is point from that cell. 'l' means left, 't' means top, and 'd' means diagonal.

## Benchmarks

As mentioned, the time and space complexity of this algorithm is O(NM). To give estimated times on how long the program will take for different sized genomes, I ran tests with differently-sized sequences. The results are as follows.

| Sequence 1/Sequence 2 (gene) | Length 1 | Length 2 | Runtime (s)
| Human/Drosophila (beta-tubulin) | ~2500 | ~1,800 | 6 |
| Prevotella/Geomicrobium (whole genome) | ~50,000 | ~50,000 |  |
| Arthroderma gypseum/Candida albican (whole genome) | ~1m | ~1m |  |

## Acknowledgments

Paper: Needleman, Saul B. & Wunsch, Christian D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". Journal of Molecular Biology. 48 (3): 443–53. doi:10.1016/0022-2836(70)90057-4. PMID 5420325
