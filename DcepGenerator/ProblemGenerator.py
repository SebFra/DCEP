#!/usr/bin/env python3
# coding: utf8

import argparse
import random
import networkx as nx
import os
import DcepGenerator.convertisseur as cv
import DcepGenerator.GraphMaker as gm


def buildSolution(copyCount):
    pathSolution = []
    for u in copyCount:
        for i in range(1, copyCount[u] + 1):
            orientation = random.choice(['R', 'F'])
            pathSolution.append(u + '_' + orientation)
    numberOfShuffle = random.choice(range(2, 11, 1))
    for i in range(1, numberOfShuffle, 1):
        random.shuffle(pathSolution)
    return pathSolution


def buildOverlapsFromSolution(pathSolution, circular):
    overlaps = list(zip(pathSolution, pathSolution[1:]))
    if circular:
        overlaps.append((pathSolution[-1], pathSolution[0]))
    return overlaps


def MarkFlowForDistanceComputation(kmer_size, pathSolution):
    flow_max = sum([int(x.split('__')[2].split('_')[0]) for x in pathSolution])
    i = pathSolution[0]
    position_dict = {}
    start_position = 0
    indexStart = 0
    for i in pathSolution:
        end_position = start_position + int(i.split('__')[2].split('_')[0]) - 1
        position_dict[indexStart] = (start_position, end_position)
        start_position = end_position - kmer_size + 1
        indexStart = indexStart + 1

    return position_dict


def buildLinksFromSolution(kmerSize, insertSize, Solution, probInsert, circular, workdir):
    position_dict = MarkFlowForDistanceComputation(kmerSize, Solution)
    DicoLinks = {}
    index = 0
    graph_solution = nx.MultiDiGraph()
    graph_solution.add_nodes_from(range(0, len(Solution), 1))
    for i in range(0, len(Solution)):
        graph_solution.node[i]['unitig'] = Solution[i]
        graph_solution.node[i]['size'] = int(Solution[i].split('__')[2].split('_')[0])

    for i in Solution:
        neighbourhood = set([j for j in range(index, len(Solution), 1) if
                             j != index and insertSize - int(Solution[j].split('__')[2].split('_')[0]) - int(
                                 Solution[index].split('__')[2].split('_')[0]) <= position_dict[j][0] -
                             position_dict[index][1] <= insertSize])
        for u in neighbourhood:
            randomNumber = random.randint(0, 100)
            if randomNumber / 100 <= probInsert:
                distance = position_dict[u][0] - position_dict[index][1]
                if ((Solution[index], Solution[u]) not in DicoLinks):
                    DicoLinks[(Solution[index], Solution[u])] = [distance]
                else:
                    DicoLinks[(Solution[index], Solution[u])].append(distance)
                graph_solution.add_edge(index, u, Type='Links', Distance=distance)
        index = index + 1

    solutionIndex = range(0, len(Solution), 1)
    overlaps = list(zip(solutionIndex, solutionIndex[1:]))
    print(overlaps)
    print(len(Solution))
    if circular:
        overlaps.append((solutionIndex[-1], solutionIndex[0]))
    for (u, v) in overlaps:
        graph_solution.add_edge(u, v, Type='Overlaps', Distance=-kmerSize + 1)
    return DicoLinks, graph_solution, position_dict


def generateLinksForStep(Links_Solution, kmerSize):
    ListLinks = []
    for (i, j) in Links_Solution:
        listDistances = Links_Solution[(i, j)]
        minDistance = min(listDistances)
        maxDistance = max(listDistances)
        maxDistance = maxDistance + 1.5 * kmerSize
        minDistance = minDistance - kmerSize
        if minDistance <= 0:
            minDistance = 0
        elif maxDistance <= 0:
            maxDistance = 0
        ListLinks.append((i, j, (int(minDistance), int(maxDistance))))
    return ListLinks


def generateOverlapsForStep(Solution, position_dict, kmerSize, circular):
    overlaps = {}
    for i in range(0, len(Solution)):
        for j in range(0, len(Solution)):
            if i != j:
                distance = position_dict[j][0] - position_dict[i][1]
                if -kmerSize < distance < int(-kmerSize / 2):
                    unitig_i = Solution[i]
                    unitig_j = Solution[j]
                    if (unitig_i, unitig_j) not in overlaps or distance <= overlaps[unitig_i, unitig_j]:
                        overlaps[unitig_i, unitig_j] = distance
    if circular:
        overlaps[((Solution[-1], Solution[0]))] = -(kmerSize - 1)
    return overlaps


def writeStepFile(copyCount, Overlaps, Links, name, workdir):
    with open(workdir + '/' + name + ".step", "w") as file:
        file.write("Contigs and coverage\n")
        file.write("====================\n\n")
        for unitig in copyCount:
            file.write(unitig + "  " + str(copyCount[unitig]) + "\n")
        file.write("\nOverlap\n")
        file.write("=======\n\n")
        for (u, v) in Overlaps:
            file.write(u + "  " + v + "  " + str(-Overlaps[(u, v)]) + "\n")
        file.write("\nLinks\n")
        file.write("=======\n\n")
        for (u, v, w) in Links:
            file.write(u + "  " + v + "  " + str(w[0]) + "  " + str(w[1]) + "\n")
    return workdir + '/' + name + ".step"


def writeSolutionReport(Solution, dicoLinksSolution, name, workdir):
    with open(workdir + '/' + name + "_solution.report", "w") as file:
        file.write('Solution of ' + name + '\n')
        nb_distances = sum([len(dicoLinksSolution[(u, v)]) for (u, v) in dicoLinksSolution])
        file.write('Nb links to satisfied : ' + str(nb_distances) + "\n")
        file.write("*********************************************************\n")
        file.write("LINKS :\n")
        for (u, v) in dicoLinksSolution:
            list_distances = dicoLinksSolution[(u, v)]
            for distance in list_distances:
                file.write(u + '\t' + v + '\t' + str(distance) + '\n')
        file.write("*********************************************************\n")
        file.write("Solution :\n")
        for u in Solution:
            file.write(u + '\n')


def generationFilesStep2(fileStep1, insertSize, kmerSize, probInsert, circular, name, workdir):
    setUnitigs = set()
    copyCount = {}
    setNonRepeatedUnitigs = set()
    with open(fileStep1, "r") as file:
        for line in file:
            lineSplit = line.split()
            setUnitigs.add(lineSplit[0])
            copyCount[lineSplit[0]] = int(lineSplit[1])
            if copyCount[lineSplit[0]] == 1:
                setNonRepeatedUnitigs.add(lineSplit[0].split('__')[0])
    Solution = buildSolution(copyCount)
    Overlaps = buildOverlapsFromSolution(Solution, circular)
    Links_Solution, graph_solution, position_dict = buildLinksFromSolution(kmerSize, insertSize, Solution, probInsert,
                                                                           circular,
                                                                           workdir)
    nx.write_graphml(graph_solution, workdir + '/' + name + '_solutionLinks.graphml')
    Links = generateLinksForStep(Links_Solution, kmerSize)
    Overlaps = generateOverlapsForStep(Solution, position_dict, kmerSize, circular)
    step_filename = writeStepFile(copyCount, Overlaps, Links, name, workdir)

    graph_file_overlaps, graph_file_links = gm.create_graphs_problem_from_step(step_filename, name)
    print("GRAPHES GENERATION : DONE")
    print("***************************************************")
    cv.graph_ml_to_data(graph_file_overlaps, graph_file_links, name)
    cv.graph_ml_to_data_solution(Solution, circular, graph_file_overlaps, graph_file_links, name)
    print("DAT GENERATION : DONE")
    print("***************************************************")
    writeSolutionReport(Solution, Links_Solution, name, workdir)
    multigraph = gm.create_multigraph_problem_simulator(graph_file_overlaps, graph_file_links, workdir)

    dico_nodes = gm.make_dico_node(multigraph)

    multigraph_DCEP = gm.build_DCEP_multigraph(multigraph)
    nx.write_graphml(multigraph_DCEP, workdir + '/Graphes/Problem/' + name + '_DCEP_multigraph.graphml')
    workdir = workdir + '/Graphes/Problem'
    cv.multigraph_DCEP_to_dat(multigraph_DCEP, workdir, name, dico_nodes)
    cv.nonRepeated_to_dat(setNonRepeatedUnitigs, workdir, name)
    multigraph_DCEP_solution = cv.multigraph_DCEP_solution(Solution, multigraph_DCEP, circular)
    cv.multigraph_DCEP_to_dat(multigraph_DCEP_solution, workdir, name + '_SOLUTION_', dico_nodes)
    exit()


def restricted_probability(x):
    x = float(x)
    if x < 0.0 or x > 1.0:
        raise argparse.ArgumentTypeError("%r not in range [0.0, 1.0]" % (x,))
    return x


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Distances and solution simulator')
    parser.add_argument('-u', "--unitigs",
                        help="unitigs file", required=True, type=str)
    parser.add_argument('-n', "--name",
                        help="name for the step 2 file", required=True, type=str)
    parser.add_argument('-i', "--insert",
                        help="insert size", required=False, type=int, default=600)
    parser.add_argument('-p', "--probability_insert",
                        help="insert probability", required=False, type=restricted_probability, default=0.8)
    parser.add_argument('-kmer', "--kmerSize",
                        help="Size of the overlaps + 1", required=True, type=int)
    parser.add_argument('-c', "--circular",
                        help="If genome is circular use this parameter", action='store_true', default=False)
    args = parser.parse_args()
    workdir = os.path.dirname(os.path.abspath(str(args.unitigs)))
    generationFilesStep2(args.unitigs, args.insert, args.kmerSize, args.probability_insert,
                         args.circular, args.name, workdir)
