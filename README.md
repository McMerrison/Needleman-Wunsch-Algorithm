# Needleman-Wunsch-Algorithm

Implementation of a DNA sequence alignment algorithm. Time and Space Complexity: O(nm) where n and m are the lengths of the sequences.

## Background:
Strands of DNA are made up of molecules called nucleotides, of which there are four types: Adenosine (A), Thymine (T), Cytosine (C), annd Guanine (G). Our genome is made up of a long sequence of millions of these nucleotides with a variable ordering. Some parts are ordered the same, making up specific genes that all of us possess. Other areas, such as those that code for aspects such as appearance, are variable even among the same species. The further two species are in similarity and genetic structure, the greater the difference in their DNA sequence. 
An important part of biological research/application is comparing two DNA sequences and assessing similarity. From an academic perspective, this can be used to determine the phylogeny of some unknown DNA sequence, or the relationship between two known species' DNA. It can also be used to show how a certain DNA sequence mutated (insertion, deletion, etc) to form a new one, and where the most likely changes occured along the chain. Assuming two sequences of size N and M (where size is measured in number of nucleotides), the time complexity of a comprehensive algorithm that checks every possible combination is O(NM). With three or more sequences, the problem would be NP-Complete. The Needleman-Wunsch algorithm utilizes dynamic programming to compare two DNA sequences and line them up based on matching nucleotides.

### How It Works
The DNA strands of interest (for now, we'll work with just two) are lined up along the row and column of a matrix, leaving an empty space at the first cell of each, as follows:

(Using sample sequences ATTCCG and ATCG)

|   | o | A | T | T | C | C | G |
| - | - | - | - | - | - | - | - |
| o |   |   |   |   |   |   |   |
| A |   |   |   |   |   |   |   |
| T |   |   |   |   |   |   |   |
| C |   |   |   |   |   |   |   |
| G |   |   |   |   |   |   |   |

Values of the matrix are initialized with a 0 in the top leftmost corner, followed by -1, -2...-n along the first row and -1, -2....-m down the first column.


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

These values can be assigned arbitrarily (a match should always be better than a mismatch/indel). We can use a generic matrix like the one above or, in real-world situations, an expert biologist would assign one for us. A match indicates that the nucleotide indexes are the same. A mismatch indicates the opposite. An insertion means one sequence is longer than the other and implies a gap. A deletion means the same, except that the sequence is shorter.

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

ATTCCG
ATC-G-



### Installing

A step by step series of examples that tell you how to get a development env running

Say what the step will be

```
Give the example
```

And repeat

```
until finished
```

End with an example of getting some data out of the system or using it for a little demo

## Running the tests

Explain how to run the automated tests for this system

### Break down into end to end tests

Explain what these tests test and why

```
Give an example
```

### And coding style tests

Explain what these tests test and why

```
Give an example
```

## Deployment

Add additional notes about how to deploy this on a live system

## Built With

* [Dropwizard](http://www.dropwizard.io/1.0.2/docs/) - The web framework used
* [Maven](https://maven.apache.org/) - Dependency Management
* [ROME](https://rometools.github.io/rome/) - Used to generate RSS Feeds

## Contributing

Please read [CONTRIBUTING.md](https://gist.github.com/PurpleBooth/b24679402957c63ec426) for details on our code of conduct, and the process for submitting pull requests to us.

## Versioning

We use [SemVer](http://semver.org/) for versioning. For the versions available, see the [tags on this repository](https://github.com/your/project/tags). 

## Authors

* **Billie Thompson** - *Initial work* - [PurpleBooth](https://github.com/PurpleBooth)

See also the list of [contributors](https://github.com/your/project/contributors) who participated in this project.

## License

This project is licensed under the MIT License - see the [LICENSE.md](LICENSE.md) file for details

## Acknowledgments

* Hat tip to anyone whose code was used
* Inspiration
* etc
