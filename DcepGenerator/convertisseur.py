#!/usr/bin/env python
# coding: utf8
import os
import sys
import pandas as pd
import networkx as nx
import GraphMaker as gm


def inv_orientation(letter):
    dico_orientation = {'R': 'F', 'F': 'R'}
    return dico_orientation[letter]


def write_id_last_block(file, list_constraints_LR):
    id_last_block = max([int(block[0]) for block in list_constraints_LR])
    file.write('param id_last_block := ' + str(id_last_block) + ';')
    write_separator_in_file(file)


def write_LR(file, list_constraints_LR):
    set_id_read = set([block[1] for block in list_constraints_LR])
    list_id_read = list(set_id_read)
    list_id_read = sorted(list_id_read, key=lambda x: int(x))
    file.write("set LR :=\n")
    for id_read in list_id_read:
        file.write(id_read + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_L_in_file(file, list_constraints_LR):
    list_of_distances = list(set([(str(block[0]), block[1], block[2], block[3]) for block in list_constraints_LR]))
    file.write("set L :=\n")
    file.write("#id_block, id_read, start_unitig, end_unitig\n")
    for (p, i, u, v) in list_of_distances:
        file.write(p + '\t' + i + '\t' + u + '\t' + v + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_distances_in_file(file, list_constraints_LR):
    list_of_distances = [(str(block[0]), block[1], block[2], block[3], block[4]) for block in list_constraints_LR]
    file.write("param distances :=\n")
    file.write("#id_block, id_read, start_unitig, end_unitig, distance\n")
    for (p, i, u, v, d) in list_of_distances:
        file.write(p + '\t' + i + '\t' + u + '\t' + v + '\t' + str(d) + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_LR_block_in_file(file, list_constraints_LR):
    file.write("set LR_BLOCK :=\n")
    file.write("#first column id_read, second id_block\n")
    set_LR_block = set([(block[1], str(block[0])) for block in list_constraints_LR])
    list_LR_block = list(set_LR_block)
    list_LR_block = sorted(list_LR_block, key=lambda x: (int(x[0]), int(x[1])))
    for (i, p) in list_LR_block:
        file.write(i + '\t' + p + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_name_in_file(file, name):
    test = name.split('/')
    if len(test) > 1:
        name = test[-1]
    file.write("param Title := " + name + ';')
    write_separator_in_file(file)


def write_start_in_file(file, graph):
    assert len(graph.nodes()) > 0, "Erreur : le graphe est vide !"
    singletons = nx.isolates(graph)
    node_length_max = max([graph.node[node]['UnitigLength'] for node in graph.nodes()])
    list_biggest_unitigs = [node for node in graph.nodes() if graph.node[node]['UnitigLength'] == node_length_max and node not in singletons]
    assert len(list_biggest_unitigs) > 0, "Il n'y a pas un seul grand unitig ?"
    file.write("param start := " + '"' + list_biggest_unitigs[0] + '"' + ';')
    write_separator_in_file(file)


def write_nodes_solution_in_file(file, nodes):
    file.write("set Vsol :=\n")
    for node in nodes:
        file.write(node + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_edges_solution_in_file(file, overlaps):
    file.write("set Osol :=\n")
    for (source, sink) in overlaps:
        file.write(source + '\t' + sink)
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_nodes_in_file(file, nodes):
    file.write("set V :=\n")
    for node in nodes:
        file.write(node + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_nodes_flow_repeat_in_file(file, nodes):
    file.write("set Flow_repeat :=\n")
    for node in nodes:
        list_successor = [n for n in nodes if
                          int(n.split('_')[0]) == int(node.split('_')[0]) and int(n.split('_')[1]) > int(
                              node.split('_')[1])]
        for successor in list_successor:
            file.write(node + '\t' + successor)
            file.write("\n")
            # file.write(node + '\t' + successor_reverse)
            # file.write("\n")
    file.write(';')
    write_separator_in_file(file)


def write_bigM_in_file(file, M):
    file.write('param M := ' + str(M) + ';')
    write_separator_in_file(file)


def write_overlaps_in_file(file, overlaps):
    file.write("set O :=\n")
    for (source, sink) in overlaps:
        file.write(source + '\t' + sink)
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_overlaps_value_in_file(file, overlaps, dico_overlaps):
    file.write("param l :=\n")
    for (source, sink) in overlaps:
        file.write(source + '\t' + sink + '\t' + str(dico_overlaps[(source, sink)][0]))
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_overlaps_value_in_file_normal(file, overlaps, dico_overlaps):
    file.write("param l :=\n")
    for (source, sink) in overlaps:
        file.write(source + '\t' + sink + '\t' + str(dico_overlaps[(source, sink)]))
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_links_in_file(file, links):
    file.write("set L :=\n")
    for (source, sink) in links:
        file.write(source + '\t' + sink)
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_l_in_file(file, overlaps, graph):
    file.write('param l :=\n')
    for (source, sink) in overlaps:
        file.write(source + '\t' + sink + '\t' + str(graph[source][sink]["Distance"]))
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_w_in_file(file, nodes, graph):
    file.write("param w :=\n")
    for node in nodes:
        file.write(node + '\t' + str(graph.node[node]["UnitigLength"]))
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_b_inf_in_file(file, links, graph):
    file.write("param b_inf :=\n")
    for (source, sink) in links:
        file.write(source + '\t' + sink + '\t' + str(graph[source][sink]['Distance_min']))
        file.write("\n")
    file.write(';')
    write_separator_in_file(file)


def write_b_sup_in_file(file, links, graph):
    file.write("param b_sup :=\n")
    for (source, sink) in links:
        file.write(source + '\t' + sink + '\t' + str(graph[source][sink]['Distance_max']))
        file.write("\n")
    file.write(';')
    write_separator_in_file(file)


def write_separator_in_file(file):
    file.write('\n')


def reverse_node(node):
    rev_node = node[0:-1] + 'F' if node[-1] == 'R' else node[0:-1] + 'R'
    return rev_node


def write_reverse_in_file(file, nodes):
    file.write("set reverse :=\n")
    for node in nodes:
        file.write(node + '\t' + reverse_node(node))
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_p_in_file(file, p):
    file.write("param p := " + str(p) + ';')
    write_separator_in_file(file)


def write_start_node(file, graph):
    source_nodes = [node for node, indegree in graph.in_degree(graph.nodes()).items() if indegree == 0]
    assert len(source_nodes) == 1, "Erreur : sources multiples"
    file.write("param start_sol := " + str(source_nodes[0]) + ';')
    write_separator_in_file(file)


def write_start(file, start):
    file.write("param start := " + str(start) + ';')
    write_separator_in_file(file)


def write_end_node(file, graph):
    sink_nodes = [node for node, outdegree in graph.out_degree(graph.nodes()).items() if outdegree == 0]
    assert len(sink_nodes) == 1, "Erreur : puits multiples"
    file.write("param end_sol := " + str(sink_nodes[0]) + ';')
    write_separator_in_file(file)


def write_w_in_file_dico(file, nodes, dico_weight):
    file.write("param w :=\n")
    for node in nodes:
        file.write(node + '\t' + str(dico_weight[node][0]))
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_b_generic(file, links, dico_b, name):
    file.write("param " + name + " :=\n")
    for (source, sink) in links:
        file.write(source + '\t' + sink + '\t' + str(dico_b[(source, sink)][0]))
        file.write("\n")
    file.write(';')
    write_separator_in_file(file)


def write_l_in_file_multi(file, multigraph):
    file.write('param l :=\n')
    for (u, v, k) in multigraph.edges(keys=True):
        if multigraph.edge[u][v][k]["Type"] == 'overlaps':
            file.write(u + '\t' + v + '\t' + str(multigraph.edge[u][v][k]["Distance"]))
            file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_w_in_file_multi(file, nodes, multigraph):
    file.write("param w :=\n")
    for node in nodes:
        file.write(node + '\t' + str(multigraph.node[node]["UnitigLength"]))
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_b_inf_in_file_multi(file, g_list, multigraph):
    file.write('param b_inf :=\n')
    for (u, v, k) in multigraph.edges(keys=True):
        if multigraph.edge[u][v][k]["Type"] == 'links' and (u, v) in g_list:
            file.write(u + '\t' + v + '\t' + str(multigraph.edge[u][v][k]["Distance_min"]))
            file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_b_sup_in_file_multi(file, g_list, multigraph):
    file.write('param b_sup :=\n')
    for (u, v, k) in multigraph.edges(keys=True):
        if multigraph.edge[u][v][k]["Type"] == 'links' and (u, v) in g_list:
            file.write(u + '\t' + v + '\t' + str(multigraph.edge[u][v][k]["Distance_max"]))
            file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_links_solution_multi(file, g_list, distance_dict):
    file.write("param solution_links :=\n")
    for (u, v) in distance_dict:
        if (u, v) in g_list:
            file.write(str(u) + '\t' + str(v) + '\t' + str(distance_dict[u, v][0]))
            file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def create_dat_minimal_problem(g_list, distance_dict, multi_graph_file, name, start):
    workdir = os.path.dirname(os.path.abspath(str(multi_graph_file)))
    multigraph = nx.read_graphml(multi_graph_file)
    nodes = multigraph.nodes()
    list_overlaps = [(u, v) for (u, v, z) in multigraph.edges(keys=True) if
                     multigraph.edge[u][v][z]['Type'] == 'overlaps']
    name_for_data_file = workdir + '/' + name + '_minimal_problem.dat'
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, name)
        write_p_in_file(file, 1)
        write_start(file, start)
        write_nodes_in_file(file, nodes)
        write_reverse_in_file(file, nodes)
        write_nodes_flow_repeat_in_file(file, nodes)
        write_overlaps_in_file(file, list_overlaps)
        write_links_in_file(file, g_list)
        write_bigM_in_file(file, 100000000)
        write_l_in_file_multi(file, multigraph)
        write_w_in_file_multi(file, nodes, multigraph)
        write_b_inf_in_file_multi(file, g_list, multigraph)
        write_b_sup_in_file_multi(file, g_list, multigraph)
        write_links_solution_multi(file, g_list, distance_dict)
    print("Fichier " + name_for_data_file + " créé")


def data_for_solution_pb(title, start, list_nodes_on_path_solution, list_edges_on_path_solution,
                         list_links_that_can_be_satisfied, dico_weight, dico_overlaps,
                         dico_b_inf, dico_b_sup, workdir):
    name_for_data_file = workdir + '/' + title + '_golden_solution.dat'
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, title)
        if start not in list_nodes_on_path_solution:
            write_start(file, gm.get_reverse(start))
        else:
            write_start(file, start)
        write_nodes_in_file(file, list_nodes_on_path_solution)
        write_overlaps_in_file(file, list_edges_on_path_solution)
        write_links_in_file(file, list_links_that_can_be_satisfied)
        write_overlaps_value_in_file_normal(file, list_edges_on_path_solution, dico_overlaps)
        write_w_in_file_dico(file, list_nodes_on_path_solution, dico_weight)
        write_b_generic(file, list_links_that_can_be_satisfied, dico_b_inf, "b_inf")
        write_b_generic(file, list_links_that_can_be_satisfied, dico_b_sup, "b_sup")
    print("Fichier " + name_for_data_file + " créé")


def graph_ml_sol_to_data(graph, name):
    name_for_data_file = name + '_golden_standart.dat'
    edges = graph.edges()
    nodes = graph.nodes()
    list_overlaps = [(source, sink) for (source, sink) in edges if graph[source][sink]["Distance"] < 0]
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_start_node(file, graph)
        write_end_node(file, graph)
        write_nodes_solution_in_file(file, nodes)
        write_edges_solution_in_file(file, list_overlaps)
    print("Fichier " + name_for_data_file + " créé")


def graph_ml_fantom_to_data(graph_file_overlaps, graph_file_links, name):
    workdir = os.path.dirname(os.path.abspath(str(graph_file_overlaps)))
    name_for_data_file = workdir + '/' + name + '.dat'
    graph_overlaps = nx.read_graphml(graph_file_overlaps)
    graph_links = nx.read_graphml(graph_file_links)
    nodes = graph_overlaps.nodes()
    list_overlaps = [(source, sink) for (source, sink) in graph_overlaps.edges()]
    list_links = [(source, sink) for (source, sink) in graph_links.edges()]
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, name)
        write_p_in_file(file, 1)
        write_start_in_file(file, graph_overlaps)
        write_nodes_in_file(file, nodes)
        write_nodes_flow_repeat_in_file(file, nodes)
        write_overlaps_in_file(file, list_overlaps)
        write_links_in_file(file, list_links)
        write_bigM_in_file(file, 100000000)
        write_l_in_file(file, list_overlaps, graph_overlaps)
        write_w_in_file(file, nodes, graph_overlaps)
        write_b_inf_in_file(file, list_links, graph_links)
        write_b_sup_in_file(file, list_links, graph_links)
    print("Fichier " + name_for_data_file + " créé")


def is_a_repetition(candidat, reference):
    candidat_split = candidat.split('_')
    reference_split = reference.split('_')
    return candidat_split[0] == reference_split[0] and candidat_split[2] == reference_split[2]


def write_T_links_compatible(file, list_links, nodes):
    file.write('set links_repeats :=\n')
    cpt_bloc = 1

    for node in nodes:
        list_preds = [u for (u, v) in list_links if v == node]
        list_succs = [v for (u, v) in list_links if u == node]
        done_preds = []
        done_succs = []
        for predecessor in list_preds:
            if len([u for u in done_preds if is_a_repetition(u, predecessor)]) == 0:
                list_repeated_preds = [u for u in list_preds if is_a_repetition(u, predecessor)]

                if len(list_repeated_preds) >= 2:
                    for u in list_repeated_preds:
                        file.write(str(cpt_bloc) + '\t' + str(u) + '\t' + str(node))
                        file.write("\n")

                    done_preds.append(predecessor)
                    cpt_bloc = cpt_bloc + 1

        for successor in list_succs:
            if len([u for u in done_succs if is_a_repetition(u, successor)]) == 0:
                list_repeated_succs = [u for u in list_succs if is_a_repetition(u, successor)]

                if len(list_repeated_succs) >= 2:
                    for u in list_repeated_succs:
                        file.write(str(cpt_bloc) + '\t' + str(node) + '\t' + str(u))
                        file.write("\n")
                    done_succs.append(successor)
                    cpt_bloc = cpt_bloc + 1
    file.write(";")
    write_separator_in_file(file)
    return cpt_bloc


def write_cpt_bloc_in_file(file, cpt_bloc):
    file.write("param cpt_bloc := " + str(cpt_bloc) + ';')
    write_separator_in_file(file)


def graph_ml_to_data(graph_file_overlaps, graph_file_links, name):
    workdir = os.path.dirname(os.path.abspath(str(graph_file_overlaps)))
    name_for_data_file = workdir + '/' + name + '.dat'
    graph_overlaps = nx.read_graphml(graph_file_overlaps)
    graph_links = nx.read_graphml(graph_file_links)
    nodes = graph_overlaps.nodes()
    list_overlaps = [(source, sink) for (source, sink) in graph_overlaps.edges()]
    list_links = [(source, sink) for (source, sink) in graph_links.edges()]
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, name)
        write_p_in_file(file, 1)
        write_start_in_file(file, graph_overlaps)
        write_nodes_in_file(file, nodes)
        write_nodes_flow_repeat_in_file(file, nodes)
        write_overlaps_in_file(file, list_overlaps)
        write_links_in_file(file, list_links)
        # write_bigM_in_file(file, 100000000)
        cpt_bloc = write_T_links_compatible(file, list_links, nodes)
        if cpt_bloc >= 2:
            cpt_bloc = cpt_bloc - 1
        write_cpt_bloc_in_file(file, cpt_bloc)
        write_l_in_file(file, list_overlaps, graph_overlaps)
        write_w_in_file(file, nodes, graph_overlaps)
        write_reverse_in_file(file, nodes)
        write_b_inf_in_file(file, list_links, graph_links)
        write_b_sup_in_file(file, list_links, graph_links)
    print("Fichier " + name_for_data_file + " créé")


def write_to_dat_generic_dim1(file, set_to_write, name_of_the_set):
    file.write("set " + name_of_the_set + " :=\n")
    for to_write in set_to_write:
        file.write(to_write + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_to_dat_generic_dim2(file, set_to_write, name_of_the_set):
    file.write("set " + name_of_the_set + " :=\n")
    for (to_write1, to_write2) in set_to_write:
        file.write(to_write1 + '\t' + to_write2 + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_to_dat_generic_dim4(file, set_to_write, name_of_the_set):
    file.write("set " + name_of_the_set + " :=\n")
    for (to_write1, to_write2, to_write3, to_write4) in set_to_write:
        file.write(to_write1 + '\t' + to_write2 + '\t' + to_write3 + '\t' + to_write4 + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_distances_satisfied_in_file(file, setDistances):
    file.write("set L :=\n")
    for (source, sink) in setDistances:
        file.write(source + '\t' + sink)
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def write_distances_satisfied_values_in_file(file, dicoDistances):
    file.write("param Distances :=\n")
    for (source, sink) in dicoDistances:
        file.write(source + '\t' + sink + '\t' + str(dicoDistances[(source, sink)][0]))
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def sol_to_dat_for_contig_generation(list_unitigs_solution, list_edges_solution, setCov, setAdj_oneSided,
                                     setAdj_twoSided, setFree, setReversibleSymetric, set_links_closed,
                                     set_links_closed_occurence_decision, dico_distances_satisfied, dico_overlaps, name,
                                     workdir):
    name_for_data_file = workdir + '/' + name + '_for_contig_generation.dat'
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, name)
        write_nodes_in_file(file, list_unitigs_solution)
        write_overlaps_in_file(file, list_edges_solution)
        write_overlaps_value_in_file(file, list_edges_solution, dico_overlaps)
        write_distances_satisfied_in_file(file, set(dico_distances_satisfied.keys()))
        write_distances_satisfied_values_in_file(file, dico_distances_satisfied)
        write_to_dat_generic_dim1(file, setCov, "Cov")
        write_to_dat_generic_dim1(file, setAdj_oneSided, "Adj_oneSided")
        write_to_dat_generic_dim1(file, setAdj_twoSided, "Adj_twoSided")
        write_to_dat_generic_dim1(file, setFree, "Free")
        write_to_dat_generic_dim1(file, setReversibleSymetric, "ReversibleSymetric")
        write_to_dat_generic_dim2(file, set_links_closed, "LinksClosed")
        write_to_dat_generic_dim4(file, set_links_closed_occurence_decision, "LinksClosedToCut")
    print("Fichier " + name_for_data_file + " créé")

def simple_sol_to_dat_for_contig_generation(list_unitigs_solution, list_edges_solution, setReversibleSymetric, set_links_closed, dico_distances_satisfied, dico_overlaps, name,
                                     workdir):
    name_for_data_file = workdir + '/' + name + '_for_contig_generation.dat'
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, name)
        write_nodes_in_file(file, list_unitigs_solution)
        write_overlaps_in_file(file, list_edges_solution)
        write_overlaps_value_in_file(file, list_edges_solution, dico_overlaps)
        write_distances_satisfied_in_file(file, set(dico_distances_satisfied.keys()))
        write_distances_satisfied_values_in_file(file, dico_distances_satisfied)
        write_to_dat_generic_dim1(file, setReversibleSymetric, "ReversibleSymetric")
        write_to_dat_generic_dim2(file, set_links_closed, "LinksClosed")
    print("Fichier " + name_for_data_file + " créé")


def LR_graph_ml_to_data(graph_overlaps, graph_file_overlaps, list_constraints_LR, name):
    workdir = os.path.dirname(os.path.abspath(str(graph_file_overlaps)))

    name_for_data_file = workdir + '/' + name + '_LR' + '.dat'
    nodes = graph_overlaps.nodes()
    list_overlaps = [(source, sink) for (source, sink) in graph_overlaps.edges()]
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, name)
        write_start_in_file(file, graph_overlaps)
        write_nodes_in_file(file, nodes)
        write_reverse_in_file(file, nodes)
        write_overlaps_in_file(file, list_overlaps)
        write_id_last_block(file, list_constraints_LR)
        write_LR(file, list_constraints_LR)
        write_L_in_file(file, list_constraints_LR)
        write_LR_block_in_file(file, list_constraints_LR)
        write_distances_in_file(file, list_constraints_LR)
        write_l_in_file(file, list_overlaps, graph_overlaps)
        write_w_in_file(file, nodes, graph_overlaps)

    print("Fichier " + name_for_data_file + " créé")


def write_start_in_file_from_DCEP(file, multigraph_DCEP):
    list_overlaps = []
    for (u, v, k) in multigraph_DCEP.edges(keys=True):
        if multigraph_DCEP.edge[u][v][k]['Type'] == 'overlaps':
            list_overlaps.append((u, v, k, multigraph_DCEP.edge[u][v][k]['Distance']))
    list_overlaps.sort(key=lambda x: x[3], reverse=True)
    file.write("param start := " + '"' + list_overlaps[0][0] + '"' + ';')
    write_separator_in_file(file)


def write_start_in_file_from_multi(file, multigraph):
    assert len(multigraph.nodes()) > 0, "Erreur : le graphe est vide !"
    node_length_max = max([multigraph.node[node]['UnitigLength'] for node in multigraph.nodes()])
    list_biggest_unitigs = [node for node in multigraph.nodes() if
                            multigraph.node[node]['UnitigLength'] == node_length_max]
    assert len(list_biggest_unitigs) > 0, "Il n'y a pas un seul grand unitig ?"
    file.write("param start := " + '"' + list_biggest_unitigs[0] + '"' + ';')
    write_separator_in_file(file)


def write_edges_for_multi(file, multigraph, name_of_the_set):
    file.write("set " + name_of_the_set + " :=\n")
    for (u, v, k) in multigraph.edges(keys=True):
        file.write(u + '\t' + v + '\t' + multigraph.edge[u][v][k]['Type'] + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_edges_value_for_multi(file, multigraph, name_of_the_param):
    file.write("param " + name_of_the_param + " :=\n")
    for (u, v, k) in multigraph.edges(keys=True):
        if multigraph.edge[u][v][k]['Type'] == 'overlaps':
            file.write(u + '\t' + v + '\t' + multigraph.edge[u][v][k]['Type'] + '\t' + str(multigraph.edge[u][v][k][
                                                                                               'Distance']) + '\n')
        elif multigraph.edge[u][v][k]['Type'] == 'links':
            file.write(u + '\t' + v + '\t' + multigraph.edge[u][v][k]['Type'] + '\t' + str(multigraph.edge[u][v][k][
                                                                                               'Distance_max']) + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_edges_value_for_multi_iteration(file, multigraph, name_of_the_param, iteration_number):
    file.write("param " + name_of_the_param + " :=\n")
    for (u, v, k) in multigraph.edges(keys=True):
        if multigraph.edge[u][v][k]['Type'] == 'overlaps':
            file.write(u + '\t' + v + '\t' + multigraph.edge[u][v][k]['Type'] + '\t' + str(multigraph.edge[u][v][k][
                                                                                               'Distance']) + '\n')
        elif multigraph.edge[u][v][k]['Type'] == 'links':
            if int(u.split('_')[1]) <= iteration_number and int(v.split('_')[1]) <= iteration_number:
                file.write(u + '\t' + v + '\t' + multigraph.edge[u][v][k]['Type'] + '\t' + str(multigraph.edge[u][v][k][
                                                                                               'Distance_max']) + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_nodes_value_in_file(file, multigraph, name_param):
    file.write("param " + name_param + " :=\n")
    for u in multigraph.nodes():
        file.write(u + '\t' + str(multigraph.node[u]['UnitigLength']) + '\n')
    file.write(";")
    write_separator_in_file(file)


def write_dico_nodes_in_file(file, dico_nodes, name_param):
    file.write("param " + name_param + " :=\n")
    for u in dico_nodes:
        file.write(u + '\t' + str(dico_nodes[u]) + '\n')
    file.write(";")
    write_separator_in_file(file)


def multigraph_DCEP_to_dat(multigraph_DCEP, workdir, name, dico_nodes):
    name_for_data_file = workdir + '/' + name + '_DCEP_multigraph' + '.dat'
    nodes = multigraph_DCEP.nodes()
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, name)
        write_start_in_file_from_DCEP(file, multigraph_DCEP)
        write_nodes_in_file(file, nodes)
        write_dico_nodes_in_file(file, dico_nodes, 'w')
        write_reverse_in_file(file, nodes)
        write_edges_for_multi(file, multigraph_DCEP, "Edges")
        write_edges_value_for_multi(file, multigraph_DCEP, "l")
    print("Fichier " + name_for_data_file + " créé")


def multigraph_to_dat(multigraph, workdir, name):
    name_for_data_file = workdir + '/' + name + '_multigraph' + '.dat'
    nodes = multigraph.nodes()
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, name)
        write_start_in_file_from_multi(file, multigraph)
        write_nodes_in_file(file, nodes)
        write_nodes_value_in_file(file, multigraph, 'w')
        write_reverse_in_file(file, nodes)
        write_edges_for_multi(file, multigraph, "Edges")
        write_edges_value_for_multi(file, multigraph, "l")
    print("Fichier " + name_for_data_file + " créé")


def write_edges_for_multi_iteration(file, multigraph, name_of_the_set, iteration_number):
    file.write("set " + name_of_the_set + " :=\n")
    for (u, v, k) in multigraph.edges(keys=True):
        if multigraph.edge[u][v][k]['Type'] != "links" or (
                int(u.split('_')[1]) <= iteration_number and int(v.split('_')[1]) <= iteration_number):
            file.write(u + '\t' + v + '\t' + multigraph.edge[u][v][k]['Type'] + '\n')
    file.write(";")
    write_separator_in_file(file)


def multigraph_to_dat_iteration(multigraph, workdir, name, iteration_number):
    name_for_data_file = workdir + '/' + name + '_multigraph_iteration_' +str(iteration_number)+ '.dat'
    nodes = multigraph.nodes()
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, name)
        write_start_in_file_from_multi(file, multigraph)
        write_nodes_in_file(file, nodes)
        write_nodes_value_in_file(file, multigraph, 'w')
        write_reverse_in_file(file, nodes)
        write_edges_for_multi_iteration(file, multigraph, "Edges", iteration_number)
        write_edges_value_for_multi_iteration(file, multigraph, "l", iteration_number)
    print("Fichier " + name_for_data_file + " créé")


def tag_solution(graph_problem, dataframeSolution):
    for i, row in dataframeSolution.iterrows():
        source = row['source']
        sink = row['sink']
        graph_problem.edge[source][sink]['solution'] = True;


def build_digraph_from_data_frame(dataframe):
    nodes = set()
    edges = set()
    graph = nx.DiGraph()

    for i, row in dataframe.iterrows():
        source = row['source']
        sink = row['sink']
        nodes.add(source)
        nodes.add(sink)
        edges.add((source, sink))
    graph.add_nodes_from(list(nodes))
    graph.add_edges_from(list(edges))
    return graph


def graph_to_graphml(graph, workdir, name):
    nx.write_graphml(graph, workdir + "/" + name + ".graphml")
    print("Le fichier " + workdir + "/" + name + ".graphml" + " a été écrit")


def csv_to_solution(graph_problem, df_solution_nodes, df_solution_edges, df_solution_holes):
    graph_solution = nx.DiGraph()
    first = True;
    circular = False
    start_save = ''
    for i, row in df_solution_nodes.iterrows():
        if first:
            start_save = row['unitig']
            first = False
        if row['source'] + row['intermediaire'] + row['sink'] >= 0.95 and row['unitig'] != "target_artificial":
            name = row['unitig']
            graph_solution.add_node(name, {'BigUnitig': graph_problem.node[name]['BigUnitig'],
                                           'UnitigLength': graph_problem.node[name]['UnitigLength']})

    for i, row in df_solution_edges.iterrows():
        if row['used'] == 1:
            source = row['source']
            if source == "target_artificial":
                source = start_save
                circular = True
            sink = row['sink']
            if sink == "target_artificial":
                sink = start_save
                circular = True
            graph_solution.add_edge(source, sink,
                                    {'Distance': graph_problem[source][sink]['Distance'], 'Flow': int(row['flow']),
                                     'Type': 'overlaps'})

    for i, row in df_solution_holes.iterrows():
        if row['satisfied'] == 1:
            source = row['source']
            if source == "target_artificial":
                source = start_save
                circular = True
            sink = row['sink']
            if sink == "target_artificial":
                sink = start_save
                circular = True
            graph_solution.add_edge(source, sink, {'Distance_min': graph_problem[source][sink]['Distance_min'],
                                                   'Distance_max': graph_problem[source][sink]['Distance_max'],
                                                   'Type': 'mate-pairs'})

    return graph_solution, circular, start_save


def workflow_after_ampl(graph_problem, df_solution_edges, df_solution_nodes, df_solution_holes, workdir, name):
    graph_solution, circular, start_save = csv_to_solution(graph_problem, df_solution_nodes, df_solution_edges,
                                                           df_solution_holes)
    if circular:
        graph_to_graphml(graph_solution, workdir, name + '_circular')
        graph_solution.remove_edges_from(
            [(source, sink) for (source, sink) in graph_solution.edges() if sink == start_save])
        graph_to_graphml(graph_solution, workdir, name + '_linear')
    else:
        graph_to_graphml(graph_solution, workdir, name)


def convert_sol_to_graphml(problemGraphml, solutionCsvNodes, solutionCsvEdges, solutionCsvHoles, name_type):

    name = problemGraphml.split('.')[0].split('/')[-1]
    name = name + '_' + name_type
    workdir = os.path.dirname(os.path.abspath(problemGraphml))



    df_solution_edges = pd.read_csv(os.path.abspath(solutionCsvEdges))
    df_solution_nodes = pd.read_csv(os.path.abspath(solutionCsvNodes))
    df_solution_holes = pd.read_csv(os.path.abspath(solutionCsvHoles))


    graph_problem = nx.read_graphml(problemGraphml)

    workflow_after_ampl(graph_problem, df_solution_edges, df_solution_nodes, df_solution_holes, workdir, name)

def writeTitle(Title, file):
    file.write("param Title :=" + str(Title))
    file.write(";")
    write_separator_in_file(file)


def writep(p, file):
    file.write("param p :=" + str(p))
    file.write(";")
    write_separator_in_file(file)


def writeStart(start, file):
    file.write("param start :=" + str(start))
    file.write(";")
    write_separator_in_file(file)


def writeM(M, file):
    file.write("param M :=" + str(M))
    file.write(";")
    write_separator_in_file(file)


def writeV(V, file):
    write_to_dat_generic_dim1(file, V, 'V')


def writeReverse(reverse, file):
    write_to_dat_generic_dim2(file, reverse, 'reverse')


def writeFlowRepeat(Flow_repeat, file):
    write_to_dat_generic_dim2(file, Flow_repeat, 'Flow_repeat')


def writeO(O, file):
    write_to_dat_generic_dim2(file, O, 'O')


def writeL(L, file):
    write_to_dat_generic_dim2(file, L, 'L')


def writel(l, file):
    file.write("param l :=\n")
    for (u, v) in l:
        file.write(u + '\t' + v + '\t' + str(int(l[(u, v)][0])))
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def writew(w, file):
    file.write("param w :=\n")
    for u in w:
        file.write(u + '\t' + str(int(w[u][0])))
        file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def writeb_inf(L, b_inf, file):
    file.write("param b_inf :=\n")
    for (u, v) in b_inf:
        if(u,v) in L:
            file.write(u + '\t' + v + '\t' + str(int(b_inf[(u, v)][0])))
            file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def writeb_sup(L, b_sup, file):
    file.write("param b_sup :=\n")
    for (u, v) in b_sup:
        if (u, v) in L:
            file.write(u + '\t' + v + '\t' + str(int(b_sup[(u, v)][0])))
            file.write("\n")
    file.write(";")
    write_separator_in_file(file)


def nb_occ(u):
    return int(u.split('_')[1])


def philter(L, nb_occ_max):
    return [(u, v) for (u, v) in L if nb_occ(u) <= nb_occ_max and nb_occ(v) <= nb_occ_max]


def create_data_for_distances_iteration(Title, p, start, V, reverse, Flow_repeat, O, L, l, w, b_inf, b_sup, workdir):
    name_for_data_file_big = workdir + '/' + Title + '_iteration0_big_SR' + '.dat'
    name_for_data_file_small = workdir + '/' + Title + '_iteration0_small_SR' + '.dat'
    L = philter(L, 0)
    V = [u[0] for u in V]
    L_small = [(u, v) for (u, v) in L if b_sup[(u, v)][0] <= 1500]
    L_big = [(u, v) for (u, v) in L if b_sup[(u, v)][0] > 1500]
    with open(name_for_data_file_small, "w") as file:
        file.write("data;\n")
        writeTitle(Title, file)
        writep(p, file)
        writeStart(start, file)
        writeV(V, file)
        writeReverse(reverse, file)
        writeFlowRepeat(Flow_repeat, file)
        writeO(O, file)
        writeL(L_small, file)
        writel(l, file)
        writew(w, file)
        writeb_inf(L_small, b_inf, file)
        writeb_sup(L_small, b_sup, file)
    with open(name_for_data_file_big, "w") as file:
        file.write("data;\n")
        writeTitle(Title, file)
        writep(p, file)
        writeStart(start, file)
        writeV(V, file)
        writeReverse(reverse, file)
        writeFlowRepeat(Flow_repeat, file)
        writeO(O, file)
        writeL(L_big, file)
        writel(l, file)
        writew(w, file)
        writeb_inf(L_big, b_inf, file)
        writeb_sup(L_big, b_sup, file)

def graph_ml_to_data_solution(Solution, circular, graph_file_overlaps, graph_file_links, name):
    workdir = os.path.dirname(os.path.abspath(str(graph_file_overlaps)))
    unitigsNumber = set([int(u.split('__')[0]) for u in Solution])
    nbOcc = {}
    SolutionOcc = []
    for unitigNumber in unitigsNumber:
        nbOcc[unitigNumber] = 0
    for u in Solution:
        u_number = int(u.split('__')[0])
        u_orientation = u.split('__')[2].split('_')[1]
        u_occ = str(u_number)+'_'+str(nbOcc[u_number])+'_'+u_orientation
        nbOcc[u_number] = nbOcc[u_number] + 1
        SolutionOcc.append(u_occ)

    overlaps_solution = list(zip(SolutionOcc, SolutionOcc[1:]))
    if circular:
        overlaps_solution.append((SolutionOcc[-1], SolutionOcc[0]))

    name_for_data_file = workdir + '/' + name + 'pathSolutionFixed.dat'
    graph_overlaps = nx.read_graphml(graph_file_overlaps)
    edges_to_remove = [(u,v) for (u,v) in graph_overlaps.edges() if (u,v) not in overlaps_solution]
    graph_overlaps.remove_edges_from(edges_to_remove)
    graph_links = nx.read_graphml(graph_file_links)
    nodes = graph_overlaps.nodes()
    list_overlaps = [(source, sink) for (source, sink) in graph_overlaps.edges()]
    list_links = [(source, sink) for (source, sink) in graph_links.edges()]
    with open(name_for_data_file, "w") as file:
        file.write("data;\n")
        write_name_in_file(file, name)
        write_p_in_file(file, 1)
        write_start_in_file(file, graph_overlaps)
        write_nodes_in_file(file, nodes)
        write_nodes_flow_repeat_in_file(file, nodes)
        write_overlaps_in_file(file, list_overlaps)
        write_links_in_file(file, list_links)
        # write_bigM_in_file(file, 100000000)
        cpt_bloc = write_T_links_compatible(file, list_links, nodes)
        if cpt_bloc >= 2:
            cpt_bloc = cpt_bloc - 1
        write_cpt_bloc_in_file(file, cpt_bloc)
        write_l_in_file(file, list_overlaps, graph_overlaps)
        write_w_in_file(file, nodes, graph_overlaps)
        write_reverse_in_file(file, nodes)
        write_b_inf_in_file(file, list_links, graph_links)
        write_b_sup_in_file(file, list_links, graph_links)
    print("Fichier " + name_for_data_file + " créé")

if __name__ == '__main__':
    problemGraphml = str(sys.argv[1])
    solutionCsvNodes = str(sys.argv[2])
    solutionCsvEdges = str(sys.argv[3])
    solutionCsvHoles = str(sys.argv[4])
    name_type = str(sys.argv[5])

    name = problemGraphml.split('.')[0]
    name = name + '_' + name_type
    workdir = os.path.dirname(os.path.abspath(problemGraphml))

    df_solution_edges = pd.read_csv(os.path.abspath(solutionCsvEdges))
    df_solution_nodes = pd.read_csv(os.path.abspath(solutionCsvNodes))
    df_solution_holes = pd.read_csv(os.path.abspath(solutionCsvHoles))

    graph_problem = nx.read_graphml(sys.argv[1])

    workflow_after_ampl(graph_problem, df_solution_edges, df_solution_nodes, df_solution_holes, workdir, name)
