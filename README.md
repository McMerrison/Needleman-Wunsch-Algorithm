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

Example RNA: AUGCGACCUACGAU <br />
Example Protein: KTVVILTNDYPQYKGMFMDINSLGDFPLVEKVHGTEAGE

To test RNA sequences (similar to DNA, except Thymine is replaced with Uracil, U instead of T), examine and try using the exRNA text files as the two input sequences. To test protein sequences (made up of various amino acids each signified by a letter, rather than the 4 nucleotides), use the exProtein text files.


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


In the cwd, run: 

```
$ perl NWA.pl test1.txt test2.txt 
```
(Note: another option is to directly populate [arg1] and [arg2] with your own text file names).
Here is the output when run on the simple example alignment we just ran through.
```
Dynamic Programming Matrix:
0   -1   -2   -3   -4   -5   -6
-1   1   0   -1   -2   -3   -4
-2   0   2   1   0   -1   -2
-3   -1   1   1   2   1   0
-4   -2   0   0   1   1   2
Dynamic Programming Matrix:
-   -   -   -   -   -   -
-   d   l   l   l   l   l
-   t   d   l   l   l   l
-   t   t   d   d   l   l
-   t   t   t   t   d   d
Alignment:

ATTCCG
AT-C-G

```

And here is an alignment for the beta-tubulin gene in homo sapiens (humans) and drosophila melanogaster (fruit flies). Though there are differences, the core DNA remains largely the same. 
```
Alignment:

CCCCAAGCCCTG-GTCAAGACGCAG--GAATG-GG-GAAGGAGCTGC-TG-GATATGGCACACACCTTAACACAAGCAGGTTAAGTACTCACTTTCCTTTGTGGTTTCACAAATGAAACCAGGATCATTTCCAATATGAAGGAAG-AGAAT-CTCTGTCAC-TGGCAAT--ATCACAGAGC-AGACTGATGGAGATT-GTTTGTGGGTATTCTATAGATTTGTGGGATTGTTGTTGGGGTAAAAAATAT-GCATG-TT--A-GAT-ACTCAA-TC-TCGA-CCTACCAAGGGCC--CCTTT-CTAGAATATCCATGATTTTTTTTGGAC-C-A-GTATCACAAAGTTCTATTTT-GATAAAACATTAACTTTTAGAAAAACAAGTAGGCTGACTTTTTCCTATTTTT-CTACACAGGTAGGAAATATGT-GCCCCGAGCA-GTCT-TGGTGGACCTAGAACCTGGGACCAT-GGACA-GCATTCGATC---TAGCAAATTAGGAGCTCTCTTTCAACCCGACAGTTTTG-TCCAT--GGTATGTTTTTCCAGAAGGTTCCACCAGGAGGAGGGGGGGATGCTTTACTGGTGCCCTTCTCTTTTCACCTTTCTTCCCCTGCTGGTTTCT-CTTTTTGGCCACAGGTA-ACTCTGGG-GCTGGCAACAACTGGGCCAAAGGCCACTACACGGAGGGAGCCGAGCTGATCGAGAATGTCCTAGAGGTGGTGAGGCACGAGAGTGAGAG--CTGTGACTGCCTGCAGGGCTTCCAGATCGTC-CACTCCCTGGGCGGGCG-CACAGGCTCCGGGATGGGCACTCTGCTCATGAAC-AAGATTAGAGAGGAGTACCCGGACCGGATCATGAATTCCT--TCAGCGTCATGCCTTCTCCCAAGGTGTCGGACACGGTGGTGGAGCCCTACAACGCGGTTCTGTCTATCCACCAGCTGATTGAGAATGCAGATGC--CTGTTTCTGCATTGACAATGAGGCCCTCTATGACATCTGCTTCCGTACCCTGAAGCTGACGACA-CCCACCTATGGGGATCTC-AACCACCTAGTGTC-CTTGACCATGAGCGGCATAACCACCTCCCTCCGGTTCCCGGGTCAGCTCAACGCAGACCTGCGCAAGCTGGCGGTGAACATGGTCCCCTTCCCCCGCCTGCACTTCTTTATGCCCGCTTT-GCCCCACTCACGG-GC-CAGGG-CAGCCAGCAGTACCGAGCCCTCTCCGTGGCCGAGCTCACCCAGCAGATGTTCGATGCCC-GCAATACCATGGCTGCCTGTGA-CCTCCGCCGTGGCCGCTACCTCACAGTGGCCTGCATT-TTCCGGGGCAAG-ATGTCCACCAAGGAAGTGGACCAGCA-ACTGCTCTCCGTGCAGACCAGGAACAGCAGCTGCTTTGTGGAGTGGATTCCC-AACAACGTCAAGGTGGCTGTCTGCGACATCCCGCCCCGG-GGGC-TGA-GCATGGCCGCCACCTTCATTGGCAACAACACGGGCATCCAAGAGATCT-TTAATAG-GG-TCTCTGAGCATTTCTCAGCCATGTTCAAAAG-GAAAG-CTTTTGTGCACTGGTACACAG-CGA----A-GG--GA--TGGA---CAT--A-----A-A-C-----GAA-----T--T-T--G-G-----G--G-A--AG----CT---G---A--A--A-----A--------T-AA-CA--------TC-----C-AT-G-ATTT---GGTAT--CC-GAG---T-A----C---CA--AC---A------AT-----T-T----TA----------A-G----A-----TG-C--C-AA--A--GCAGT-TC-TA-GA-----G--G----A-AGAT-----G--A-AG-AG--G--------------------TCA-C------GGA--GGA-GG-CA-G-A-A--AT--GGA-GC-CA-G---A---AGA---T-AAG-GGACAT-T----AA-C---TG-T---G--AG-----AG-A---AGCTG-T-----G-C--C-------GC-------G--G-AG--TCGC---T-T-ACA---GA-A-CA-G--T-TTC-TCATTA-G---AT----G--AGT-GT----T--TCT-C-CTG----C---AG-C-A--CT---C------C-A-A-A-AC-C----CAC---------TCT---GCA---CT-GCAGC-A----CA---G-TG---AATG-AT--A----------TG-CACTC-A--C----C---A--TTA------G-C-TTCG-A--C-AC-AGG------G-A-C-T--G-AG--GGAG-ACA-G------GT---G----GG-----GA------GCAG----CTGA---C-AG---G-C-A---T-TAGG--GTCT----T---G-CT-GAC---A-T--C-T--A---CTA-ACCTTGAA-GAGTT--TG-AT-G-T--TC-AGTG-CATACT-T-ATTA-A-C-T-TA-AAAA---AAT-----A---GCA-----AATT-TATTGTA-A-AGT----GG--A---TC----C-CTTTGT-TTCA-----A--A--------GT-GTTTGCCA-GG-------CAT-----CC-AG-AC-TA-C-AGTG-T---GG-A----TTT--G-CAG---------GGAGCCCACT--------
TC--A-GT--TGTGTCTTGT-GC-GTTG--TGCGGTGAAGTTGC-GCCTTCGTTCTG--A-A-ACCT-AATA-A-GCTG---A-GAA-TCAGAA-CCTCTG-G-----ACATCTGAA-C--G--TC-TGA--AA-ATCCAGCA-GCAGC-TACTCGGAC-CGTG-CGGTGCATC-CTTCTCCAGACTGA-G-AC-TAAGATCC-G---ATTCCA--GCCTT-----A--GCC-TTCGACCCCACAA-ATCGCC-GCTTTCAAGATGAGAGAAATCGT-GAACCTGC-A-GG-CCGGCCAGTGCG-GCA-A-CCAA-ATC------GG-CGCTAAGT-TCTGGGAGATC-ATTTCCGAGGAG-CACGG-CATC--GACAG-CAA-T-GGC--A-T----C-TACGTGGGCGACA--G-T-G-A--TCTGCAGCTG-GAGCGCGTCAGTG-TCTAC-TACAAC--GA-AGCATCGGTCACGCGGTCG-TCGGGTGGCAAGT-ACGTGC-C-C----AGGGC--CA-TCCTGCTCGATCTGG-A-GC----CCGGAA----CCAT--GGAGTCGG------TGCGTTCC-GGT-CCGTACGGA---CAACT--CTTCCG--GCCGGA--CAACTTCGTGT--AC-GG-ACAGTC-GGGAGCGGGCAACAACTGGGCCAAGGGTCACTACACCGAGGGCGCCGAGCTGGTGGACAATGTCCTGGACGTGGTCCG-CAAG-GAGTGCGAGAACTGCGACTGCCTGCAGGGCTTCCA-ATTGACGCACTCGCTGGGCGG-CGGCACTGGGTCCGGAATGGGCACCCTGCTGATCT-CGAAGATCCGCGAGGAGTACCCCGACCGCATCATGAACACCTACTC-G-GTGGTGCCATCGCCCAAGGTGTCGGACACCGTGGTGGAGCCCTACAACGCCACCCTGTCGATCCACCAGCTGGTGGAGAACACAGACGAGAC-GTA-CTGCATCGACAACGAGGCGCTGTACGACATCTGCTTCCGGACGCTGAAGGTGTCGA-ATCCCAGCTACGGAGA-CTTGAACCACCTGGTCTCGCT-GACCATGTCCGGGGTGACCACCTGCCTGCGTTTCCCCGGCCAGCTGAACGCCGATCTGCGCAAGCTGGCGGTCAACATGGTTCCATTCCCGCGTCTCCACTTCTTCATGCCCGGATTCGCGCCGCTCACCTCGCGC-GGATC-GC-AGCAGTACCGCGCCCTCACCGTTCCCGAACTGACCCAGCAGATGTTCGACGCCAAG-AACATGATGGCCGCCTGTGATCCAC-GCCACGGTCGCTACCTCACGGTGGCC-GCCGTCTTCCGCGGCC-GCATGTCCATGAAGGAGGTGGACGAGCAGA-TGCTGGCGGTGCAGAACAAGAACAGCTCCTACTTCGTGGAGTGGAT-CCCGAACAATGTGAAGACCGCGGTGTGCGACATCCCGCC--GAAGGGCCTGAAG-ATGTCCTCAACGTTCATCGGCAACACCACGGCCATCCAGGAG--CTGTTCA-AGCGGATCTCCGAGCAGTTCTCGGCCATGTTCCG--GCGCAAGGCCTTCCTGCACTGGTACACCGGCGAGGGCATGGACGAGATGGAGTTCACCGAGGCGGAGAGCAACATGAACGACCTGGTGTCCGAGTACCAGCAGTACCAGGAGGCCACCGCCGACGACGAGTTCGACCCGGAGGTCAATCAGGAGGAGGTCGAGGGCGATTGTATTTAATGGTGTGGCCAGAGGAATGAGGGGCGTGCGTGACGGCAGGGGGGATAGGGATGTCCGCTTTCTGCGCGGCATGCGCCACCAGCTGGCGGCCAAGCACTGCAGCATCCTGTGACCCCCGCCGTTGGACACACACACAGTCACAGCAGACGCAGACCACACACCACCACAATCATCACACACGGACAGGAAGGACATGCAGACGAACAGGAAGGTCAAGCGTAGTTAGAACATCAAACGGAAATGTTCCCAAGCCCATGCTTCCGCTACTCCCCATCATCCAGCTCATCCATCGTCATCCATCATCGCCCATCCCGCCGCATCATCTCACATCTCACATCTGACATCTCGCATCTTCATCTTTTTGTTCATTTTCGGCACTTGTCACCTCCTCTGCGCTATATACTTTATACTATACTATACTATATACCATATACACACTTTACACGGAGGCACGTCTATCGCATTTCTCGCAGCCATTTTCACCCGCTTTCCAACCCACCCACCACCACCCATGTCACTCCATTCATAACCCCAAGTTATCCCACGTCATTCTTATACTACTATTTATATCGCAACGTCAGCAGAAGGACCACAAGTATCTAGTATTGCAAAGGATATCGAATTTATGCAGAAATCTCAAAACGAGAATGTCCAGCGTCTAGGAAGTCGGCACTCAAGACTTGGCGGCACTTACGTTTAGGGCTATACGATAAACGAGCTCATCCAAAGATAATCGAACGACA-ACTCTCATGAGATCCTCTATAAAAGGCAATTTCTGAAAAGCATCTGCAATTGTATT-TACATAATCTTAGGGCACAATCTATTCACTTTGTATTCGTGTGCAGTACTTCTCTTGTAGTTTTATATGGATAAAAACATTTGTACCCAATAAATAACTAATAATAATGGCAAACATTTTTGACAAATAATAAATGGAATTTATTTAAAATTT

Score: 341
Runtime: 6 seconds

```

When the run is complete, the program will print the original sequence, the scoring matrix (for shorter sequences), the traceback matrix (for shorter sequences), and the alignment along with final score (positive usually means good).
 (In the traceback matrix, each letter in the cell tells you where the arrow is point from that cell. 'l' means left, 't' means top, and 'd' means diagonal).

## Benchmarks

As mentioned, the time and space complexity of this algorithm is O(NM). To give estimated times on how long the program will take for different sized genomes, I ran tests with differently-sized sequences. The results are as follows.

| Sequence 1/Sequence 2 (gene) | Length 1 | Length 2 | Runtime (s)
| - | - | - | - |
| Human/Drosophila (beta-tubulin) | ~2500 | ~1,800 | 6 |
| Prevotella/Geomicrobium (whole genome) | ~50,000 | ~50,000 |  |
| Arthroderma gypseum/Candida albican (whole genome) | ~1m | ~1m |  |

## Other Methods

Various other alignment algorithms exist based on the nature of the comparison (interspecies vs. ancestor), the number of sequences, the length of sequences, and various other parameters. Most are a variation or improvement on the dynamic programming approach employed by the Needleman-Wunsch Algorithm.

### Global vs. Local

THE NW algorithm is meant for global alignment, meaning we wish to align every possible nucleotide in the entire sequence. This is useful when comparing similar species or the same gene from different species. Local Alignment is simply a variation that only seeks to find regions of similarities when comparing sequences. This is useful for finding shared genes when scanning entire genomes. The method is called the Smith-Waterman Algorithm. It is essentially the same as NW, except during the traceback step we start at the highest score value on the matrix, rather than always at the bottom rightmost cell.

### Pairwise Alignment

The simplest approach to pairwise alignment (global or local alignment of 2 sequences) is a dot-matrix, in which we design a matrix similar to the NW algorithm, except rather than populate cells with scores, we simply mark those that intersect at matching nucleotides along the top row and left column. The clarity of a diagonal marking indicates how closely related the sequences are, though there is a significant amount of wasted space, and lack of quantitative data. The scoring scheme emplyed by NW is an improvement on this rough method.

Faster word-related methods employ a heuristic that isn't guaranteed to be correct, but performs better than NW. They are often used in official gene databses like FASTA and BLAST (where every user search is essentially an pairwise alignment with multiple possible results). For each query, subsequences are identified and matched against the databse. The character indexes are compared to get an offset, and candidate matches to the full sequence. By matching small subsequences first, we can eliminate a vast majority of candidate sequences int he data base that do not contain the word (and by extension, cannot contain the query).

### Multiple Sequence Alignment

As mentioned, multiple sequence alignment is an NP-Complete problem. Traditional algorithms used for pairwise alignment can be run on multiple sequences, but the problem size and runtime will grow exponentially. Heuristic methods such as Hidden Markov Models, simulated annealing and genetic algorithms (interestingly enough, based on the very principle of chromosomal mutation, applied to a wide variety of computational problems) have also been applied. Because these more computationally efficient methods aren't always accurate, comparisons between 3 or 4 sequences are often done through traditioanl dynamic programming, since modern hardware can still perform the alignment in a reasonable amount of time, well worth an accurate alignment.

#### T-Coffee

Tree-based Consistency Objective Function for Alignment Evaluation is software that performs MSA through a series of pairwise alignments (dynamic programming), arranging sequences by the alignment score (since the same scoring scheme is used on all alignments, the actual numbers are arbitrary and relative). They are clustered and sorted by sequences that are most closely related, followed by less closely-related, and constructed into a phylogentic tree ('phylogeny' refers to evolutionary family). The time complexity is essentially O(mn^2), where n^2 refers to the pairwise alignment and m refers to the number of sequences. While not ideal, it is still polynomial and produces and accurate alignment, adaptable to numerous databases and software.

#### Iterative Methods

The main drawback of tree methods (such as T-Coffee) is the dependence on an accurate initial  alignment. The entire tree structure depends on a good initial matching to form the cluster of related species, which is then used as a relative baseline for the rest of the tree. As an alternative, iterative methods use an objective function, optimized using a scoring scheme. It starts with an initial global alignment, then realigns subsequences, which then form the basis for the next alignment.

#### Profile Analysis




## Acknowledgments

Paper: Needleman, Saul B. & Wunsch, Christian D. (1970). "A general method applicable to the search for similarities in the amino acid sequence of two proteins". Journal of Molecular Biology. 48 (3): 443–53. doi:10.1016/0022-2836(70)90057-4. PMID 5420325

Databse: https://blast.ncbi.nlm.nih.gov/Blast.cgi?PAGE_TYPE=BlastSearch

Paper: Notredame C, Higgins DG, Heringa J (2000-09-08). "T-Coffee: A novel method for fast and accurate multiple sequence alignment". J Mol Biol. 302 (1): 205–217. doi:10.1006/jmbi.2000.4042. PMID 10964570.
