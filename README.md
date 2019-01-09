## Introduction
The DCEP problem (for Distance Constraint elementary path) has been presented for the first time during the CTW 2018.
The problem consists to find a path that'll satisfy the maximum number of distances in a given Set.
Initially it's coming from the need to solve simultaneously two well known NP-Complet genomic problem called scaffolding and gap-filling.
This simulator is a tool to produce instances of this problem similar of the ones we see in genomics and to benchmark different approaches to solve them.

## Usage
There are two main steps to produce an instance of the DCEP problem :

#First we need to produce the non-oriented unitigs :
```sh
UnitigsGenerator.py -Nnr 2 -Nr 4 -MaxR 3
```

### Required
*-Nnr : number of non repeated unitig
*-Nr : number of repeated unitig

The non repeated unitig will appear only once in the solution.
The repeated unitig will appear between two and MaxR in the solution.

### Optionnal
* -MaxR   : maximum number of occurence for a repeated unitig in the solution. By default : 3
* -name   : name of the unitigs file. By default : NO_NAME
* -MinWnr : minimum length of sequence for a non repeated unitig
* -MaxWnr : maximum length of sequence for a non repeated unitig
* -MinWr  : minimum length of sequence for a repeated unitig
* -MaxWr  : maximum length of sequence for a repeated unitig

### OUTPUT
A file name_step1.txt that contains all the unitigs, their length and their copy count.
Example :
'''
1__len__5598	1
2__len__1977	1
3__len__144	    2
4__len__145	    2
5__len__160	    2
6__len__160	    2

'''

How to read it ? For the first line : "The unitig 1 has a length of 5598 characters and appeared only once in the solution".

# Last : Distances, problem and solution generation
```sh
ProblemGenerator.py -s1 NO_NAME_step1.txt -kmer 101 -n testClassicSmall

```
### Required
*-u : unitigs file (generated with UnitigsGenerator.py)
*-n : name of the problem generated
* -kmer : overlaps maximum = kmer - 1, must be superior to the smallest unitig.

The non repeated unitig will appear only once in the solution.
The repeated unitig will appear between two and MaxR in the solution.

### Optionnal
* -i : insert size
* -p : probability of insert (distance generation)
* -c : if present -> produce a circular solution otherwise linear solution.

### OUTPUT
The solution of the instance is in the name_solution.report file (use any text editor to read it).
Example :
```
Solution of testClassicSmall
Nb links to satisfied : 8
*********************************************************
LINKS :
1__len__5598_F	4__len__145_R	-41
1__len__5598_F	6__len__160_R	3
1__len__5598_F	5__len__160_R	62
5__len__160_F	1__len__5598_F	-56
4__len__145_F	1__len__5598_F	-100
2__len__1977_F	1__len__5598_F	3
3__len__144_R	1__len__5598_F	1922
3__len__144_R	1__len__5598_F	1879
1__len__5598_F	6__len__160_F	-100
*********************************************************
Solution :
3__len__144_R
3__len__144_R
2__len__1977_F
5__len__160_F
4__len__145_F
1__len__5598_F
6__len__160_F
4__len__145_R
6__len__160_R
5__len__160_R
```

A file name_solutionLinks.graphml that contains the graph of the solution. You can use cytoscape to read it.

In the directory Graphes you'll find a directory Problem with all the files regarding the problem that has been generated.
The file name.dat contains a classical representation of the problem (with the length of the unitig on the nodes and the two types of edges :
links (distances) with positive weights on and the overlaps with negative weights on.
The other version called name_DCEP_multigraph.dat is another representation of the same problem with the weight only on the edges. (and always positives).
You can use both but be carefull to use the right models to solve them.



