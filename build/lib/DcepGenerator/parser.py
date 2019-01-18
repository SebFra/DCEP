#!/usr/bin/env python

import re
import sys
import DcepGenerator.graph_tools as graph_tools


def parse_nodes(contigs_list):
    # print(nodes_str)
    list_nodes = []
    pattern = '[0-9]+'
    #  pattern = '\n(\d+)[_a-z]+(\d+)\s+\\2\s+(\d+)(?=\n)'
    for contigs in contigs_list:
        list_content_nodes = re.findall(pattern, contigs, re.MULTILINE)
        if len(list_content_nodes) == 3:
            list_nodes.append(
                graph_tools.Node(int(list_content_nodes[0]), int(list_content_nodes[1]), int(list_content_nodes[2])))
        elif len(list_content_nodes) == 4:
            list_nodes.append(
                graph_tools.Node(int(list_content_nodes[0]), int(list_content_nodes[1]), int(list_content_nodes[3])))

    return list_nodes


def parse_overlaps(overlaps_list):
    # print(edges_str)
    #  pattern = '\n(\d+)__len__(\d+)__(\w)\s+(\d+)__len__(\d+)__(\w)\s+(-?\d+)(?=\n)'
    #  pattern = '[0-9]+[A-Za-z0-9_]*'
    pattern = '[0-9]+'
    list_edge = []

    for composant_overlaps in overlaps_list:
        composants_edge = composant_overlaps.split()
        contig_source = composants_edge[0]
        part1 = re.findall(pattern, contig_source, re.MULTILINE)
        index_contig_source = int(part1[0])
        orientation_contig_source = graph_tools.INV_ORIENTATION[contig_source[-1]]
        contig_sink = composants_edge[1]
        part2 = re.findall(pattern, contig_sink, re.MULTILINE)
        index_contig_sink = int(part2[0])
        orientation_contig_sink = graph_tools.INV_ORIENTATION[contig_sink[-1]]
        distance = -abs(int(composants_edge[2]))
        list_edge.append(graph_tools.Edge(index_contig_source, index_contig_sink, distance, orientation_contig_source,
                                          orientation_contig_sink))

    return list_edge


def parse_link(link_list_str):
    # print(solution_str)
    #  pattern = '[0-9]+__len__[0-9]+[_]{1,2}[FR]'
    #  pattern = '[0-9]__len__[0-9]+_[FR][ ]*[0-9]__len__[0-9]+_[FR][ ]*[-]?[0-9]+'
    pattern = '[0-9]+'
    list_link = []

    for link_str in link_list_str:
        composants_link = link_str.split()
        print(composants_link)
        contig_source = composants_link[0]
        contig_sink = composants_link[1]
        part1 = re.findall(pattern, contig_source, re.MULTILINE)
        contig_source_index = int(part1[0])
        contig_source_orientation = graph_tools.INV_ORIENTATION[contig_source[-1]]
        part2 = re.findall(pattern, contig_sink, re.MULTILINE)
        contig_sink_index = int(part2[0])
        contig_sink_orientation = graph_tools.INV_ORIENTATION[contig_sink[-1]]
        distance_min = int(composants_link[2])
        distance_max = int(composants_link[3])
        list_link.append(
            graph_tools.Link(contig_source_index, contig_source_orientation, contig_sink_index, contig_sink_orientation,
                             distance_min, distance_max))
    return list_link


def complete_list_edges(list_edges):
    set_edges = set(list_edges)
    for edge in list_edges:
        set_edges.add(edge.get_symetrique())
    return list(set_edges)


def parse_file(filepath):
    #  pattern = '([ ]*[0-9]+[a-zA-Z0-9_ \-+\n]+?\n)\n'
    #  pattern = '^[0-9][A-Za-z0-9_]+ [ ]* [0-9][A-Za-z0-9_]+ [ ]* [0-9]+$'
    pattern_contig = '^[0-9]+__len__[0-9]+[ ]*[0-9]*[ ]*[0-9]+[ ]*[0-9]*$'
    pattern_overlap = '^[0-9]+__len__[0-9]+[_]{1,2}[FR][ ]*[0-9]+__len__[0-9]+[_]{1,2}[FR][ ]*[-]?[0-9]+$'
    pattern_link = '^[0-9]+__len__[0-9]+[_]{1,2}[FR][ ]*[0-9]+__len__[0-9]+[_]{1,2}[FR][ ]*[0-9]+[ ]+[0-9]+$'

    with open(filepath) as fi:
        text = fi.read()

    #  print(parts)
    #  for b in parts:
    #      print(b)

    contigs_part = re.findall(pattern_contig, text, re.MULTILINE)
    overlap_part = re.findall(pattern_overlap, text, re.MULTILINE)
    link_part = re.findall(pattern_link, text, re.MULTILINE)
    #  print("-------------------------------------------------------------------")
    #  print("Les contigs : ")
    #  print(contigs_part)
    #  print("-------------------------------------------------------------------")
    #  print("Les overlaps : ")
    #  print(overlap_part)
    #  for a in overlap_part:
    #      print(a)
    #  print("-------------------------------------------------------------------")
    print("Les links : ")
    print(link_part)
    print("-------------------------------------------------------------------")

    node_pool = parse_nodes(contigs_part)
    print("NODES VERIFICATION APRES CREATION DES OBJETS")
    for node in node_pool:
        print(node)
    print("-------------------------------------------------------------------")
    edges_pool = parse_overlaps(overlap_part)
    print("EDGES VERIFICATION APRES CREATION DES OBJETS")
    for edge in edges_pool:
        print(edge)
    print("-------------------------------------------------------------------")
    edges_pool_completed = complete_list_edges(edges_pool)
    print("EDGES VERIFICATION APRES COMPLETION")
    for edge in edges_pool_completed:
        print(edge)
    print("-------------------------------------------------------------------")
    links_pool = parse_link(link_part)
    links_pool_completed = complete_list_edges(links_pool)
    print("LINKS VERIFICATION APRES COMPLETION")
    for link in links_pool_completed:
        print(link)
    print("-------------------------------------------------------------------")
    print("AJOUT DES LINKS")
    for link in links_pool_completed:
        edges_pool_completed.append(link.get_edge())
    print("-------------------------------------------------------------------")
    for edge in edges_pool_completed:
        print(edge)

    return (node_pool, edges_pool_completed, links_pool_completed)

def parse_file_with_edges_division(filepath):
    #  pattern = '([ ]*[0-9]+[a-zA-Z0-9_ \-+\n]+?\n)\n'
    #  pattern = '^[0-9][A-Za-z0-9_]+ [ ]* [0-9][A-Za-z0-9_]+ [ ]* [0-9]+$'
    pattern_contig = '^[0-9]+__len__[0-9]+[ ]*[0-9]*[ ]*[0-9]+[ ]*[0-9]*$'
    pattern_overlap = '^[0-9]+__len__[0-9]+[_]{1,2}[FR][ ]*[0-9]+__len__[0-9]+[_]{1,2}[FR][ ]*[-]?[0-9]+$'
    pattern_link = '^[0-9]+__len__[0-9]+[_]{1,2}[FR][ ]*[0-9]+__len__[0-9]+[_]{1,2}[FR][ ]*[0-9]+[ ]+[0-9]+$'

    with open(filepath) as fi:
        text = fi.read()

    # print(parts)
    #  for b in parts:
    #      print(b)

    contigs_part = re.findall(pattern_contig, text, re.MULTILINE)
    overlap_part = re.findall(pattern_overlap, text, re.MULTILINE)
    link_part = re.findall(pattern_link, text, re.MULTILINE)
    #  print("-------------------------------------------------------------------")
    #  print("Les contigs : ")
    #  print(contigs_part)
    #  print("-------------------------------------------------------------------")
    #  print("Les overlaps : ")
    #  print(overlap_part)
    #  for a in overlap_part:
    #      print(a)
    #  print("-------------------------------------------------------------------")
    print("Les links : ")
    print(link_part)
    print("-------------------------------------------------------------------")

    node_pool = parse_nodes(contigs_part)
    print("NODES VERIFICATION APRES CREATION DES OBJETS")
    for node in node_pool:
        print(node)
    print("-------------------------------------------------------------------")
    overlaps_pool = parse_overlaps(overlap_part)
    print("Overlaps VERIFICATION APRES CREATION DES OBJETS")
    for edge in overlaps_pool:
        print(edge)
    print("-------------------------------------------------------------------")
    edges_overlaps_pool_completed = complete_list_edges(overlaps_pool)
    print("Overlaps VERIFICATION APRES COMPLETION")
    for edge in edges_overlaps_pool_completed:
        print(edge)
    print("-------------------------------------------------------------------")
    links_pool = parse_link(link_part)
    links_pool_completed = complete_list_edges(links_pool)
    edges_links_pool_completed = [link.get_edge() for link in links_pool_completed]
    print("LINKS VERIFICATION APRES COMPLETION")
    for link in links_pool_completed:
        print(link)
    for edge in edges_overlaps_pool_completed:
        print(edge)

    return node_pool, edges_overlaps_pool_completed, edges_links_pool_completed

if __name__ == '__main__':
    parsed_data = parse_file(sys.argv[1])

    # print('==========================================================')

#  print_reformatted_data(parsed_data)

# print(parse_file('acorus.txt'))
