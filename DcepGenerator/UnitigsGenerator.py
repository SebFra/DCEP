#!/usr/bin/env python3
# coding: utf8

import argparse
import random


def makeWdict(list_of_unitig, min_w, max_w, w):
    for u in list_of_unitig:
        w[u] = random.randrange(min_w, max_w, 1)
    return w


def makeRepDict(list_of_repeated_unitig, MaxNumberOfRepetition):
    repetition_dict = {}
    for u in list_of_repeated_unitig:
        repetition_dict[u] = random.randrange(2, MaxNumberOfRepetition, 1)
    return repetition_dict


def generationFileStep1(NumberNonRepeatedUnitigs, NumberRepeatedUnitigs, MaxNumberOfRepetition,
                        name, MinimumWeightNonRepeated,
                        MaximumWeightNonRepeated, MinimumWeightRepeated,
                        MaximumWeightRepeated):
    list_of_non_repeated_unitig = range(1, NumberNonRepeatedUnitigs + 1, 1)
    list_of_repeated_unitig = range(NumberNonRepeatedUnitigs + 1, NumberNonRepeatedUnitigs + NumberRepeatedUnitigs + 2,
                                    1)
    w = {}
    w = makeWdict(list_of_non_repeated_unitig, MinimumWeightNonRepeated, MaximumWeightNonRepeated, w)
    w = makeWdict(list_of_repeated_unitig, MinimumWeightRepeated, MaximumWeightRepeated, w)
    repetition_dict = makeRepDict(list_of_repeated_unitig, MaxNumberOfRepetition)
    with open(name+'_step1.txt', "w") as file:
        for u in list_of_non_repeated_unitig:
            file.write(str(u)+'__len__'+str(w[u])+'\t'+'1'+'\n')
        for u in list_of_repeated_unitig:
            file.write(str(u)+'__len__'+str(w[u])+'\t'+str(repetition_dict[u])+'\n')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Unitig and repetitions simulator')
    parser.add_argument('-Nnr', "--Number_Non_Repeated_Unitigs",
                        help="Number of non repeated unitigs that you want to have in the solution", required=True, type=int)
    parser.add_argument('-Nr', "--Number_Repeated_Unitigs",
                        help="Number of repeated unitigs that you want to have in the solution", required=True, type=int)
    parser.add_argument('-MaxR', "--Max_Number_Of_Repetition",
                        help="Maximum of occurence for a repeated unitig", required=False, default=3, type = int)
    parser.add_argument('-n', "--name",
                        help="Name of the instance", required=False, default="NO_NAME", type=str)
    parser.add_argument('-MinWnr', "--Minimum_Weight_Non_Repeated",
                        help="Minimum Weight for a non repeated unitig", required=False, default=200, type = int)
    parser.add_argument('-MaxWnr', "--Maximum_Weight_Non_Repeated",
                        help="Maximum Weight for a non repeated unitig", required=False, default=50000, type = int)
    parser.add_argument('-MinWr', "--Minimum_Weight_Repeated",
                        help="Minimum Weight for a repeated unitig", required=False, default=102, type = int)
    parser.add_argument('-MaxWr', "--Maximum_Weight_Repeated",
                        help="Maximum Weight for a repeated unitig", required=False, default=199, type = int)
    args = parser.parse_args()
    generationFileStep1(args.Number_Non_Repeated_Unitigs, args.Number_Repeated_Unitigs, args.Max_Number_Of_Repetition,
                        args.name,
                        args.Minimum_Weight_Non_Repeated, args.Maximum_Weight_Non_Repeated, args.Minimum_Weight_Repeated,
                        args.Maximum_Weight_Repeated)
