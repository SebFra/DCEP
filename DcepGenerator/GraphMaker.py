#!/usr/bin/env python

import networkx as nx
import itertools
import os
import parser

LIMIT_SIZE = 100


def get_reverse_from_scaffold(scaffold):
    scaffold_reversed = list(scaffold)
    scaffold_reversed.reverse()
    scaffold_reversed_complemented = []
    for unitig in scaffold_reversed:
        scaffold_reversed_complemented.append(get_reverse(unitig))
    return scaffold_reversed_complemented


def build_Digraph_problem(max_links, list_overlaps, dico_nodes_weight, List_scaffolds):
    G = nx.DiGraph()
    G.add_edges_from(list_overlaps)
    for node in G.nodes():
        G.node[node]['content'] = str(node)
        G.node[node]['weight'] = dico_nodes_weight[node][0]
        if dico_nodes_weight[node][0] > max_links:
            G.node[node]['Big'] = True
        else:
            G.node[node]['Big'] = False
    return compressToMetaDiGraph(G, List_scaffolds)


def contracted_nodes(G, u, v, self_loops=True):
    H = G.copy()
    if H.is_directed():
        in_edges = ((w, u, d) for w, x, d in G.in_edges(v, data=True)
                    if self_loops or w != u)
        out_edges = ((u, w, d) for x, w, d in G.out_edges(v, data=True)
                     if self_loops or w != u)
        new_edges = itertools.chain(in_edges, out_edges)
    else:
        new_edges = ((u, w, d) for x, w, d in G.edges(v, data=True)
                     if self_loops or w != u)
    v_name = str(v)
    H.remove_node(v)
    H.add_edges_from(new_edges)

    H.node[u]['content'] += '__' + v_name

    return H


def compressToMetaDiGraph(multiDiGraph, list_scaffolds):
    multiDiGraph_copy = multiDiGraph.copy()
    compteur_scaffold = 0
    for scaffold in list_scaffolds:
        if len(scaffold) > 1:
            first_node = scaffold[0]

            for node in scaffold[1:]:
                multiDiGraph_copy_temp = multiDiGraph_copy.copy()
                multiDiGraph_copy = contracted_nodes(multiDiGraph_copy_temp, first_node, node, self_loops=False)

            multiDiGraph_copy = nx.relabel_nodes(multiDiGraph_copy, {first_node: "S" + str(compteur_scaffold) + 'F'})
            scaffold_reversed_complemented = get_reverse_from_scaffold(scaffold)
            first_node = scaffold_reversed_complemented[0]
            for node in scaffold_reversed_complemented[1:]:
                multiDiGraph_copy_temp = multiDiGraph_copy.copy()
                multiDiGraph_copy = contracted_nodes(multiDiGraph_copy_temp, first_node, node, self_loops=False)

            multiDiGraph_copy = nx.relabel_nodes(multiDiGraph_copy, {first_node: "S" + str(compteur_scaffold) + 'R'})
        compteur_scaffold += 1
    return multiDiGraph_copy


def compressToMetaGraph(multiDiGraph, list_scaffolds):
    multiDiGraph_copy = multiDiGraph.copy()
    compteur_scaffold = 0
    for scaffold in list_scaffolds:
        if len(scaffold) > 1:
            first_node = scaffold[0]
            for node in scaffold[1:]:
                multiDiGraph_copy_temp = multiDiGraph_copy.copy()
                multiDiGraph_copy = contracted_nodes(multiDiGraph_copy_temp, first_node, node, self_loops=False)
        multiDiGraph_copy = nx.relabel_nodes(multiDiGraph_copy, {first_node: "S" + str(compteur_scaffold) + 'F'})
        compteur_scaffold += 1
    return multiDiGraph_copy


def generate_meta_graph(max_links, dico_nodes_weight, Dico_overlaps, list_links_satisfied, Path_solution,
                        List_scaffolds, workdir):
    # GRAPH GENERATION
    multi_graph = nx.MultiDiGraph()
    multi_graph.add_path(Path_solution, type="overlaps", path=True)
    for node in multi_graph.nodes():
        multi_graph.node[node]['content'] = str(node)
        multi_graph.node[node]['weight'] = dico_nodes_weight[node][0]
        multi_graph.node[node]['order'] = Path_solution.index(node) / 1.0
        if dico_nodes_weight[node][0] > max_links:
            multi_graph.node[node]['Big'] = True
        else:
            multi_graph.node[node]['Big'] = False
    # Make circular
    print(Path_solution[0], Path_solution[-1], max_links)
    multi_graph.add_edge(Path_solution[-1], Path_solution[0], type="overlaps", path=True)
    multi_graph.add_edges_from(list_links_satisfied, type="links", path=False, satisfied=True)
    graph_to_graphml(multi_graph, workdir, "solution")

    meta_graph_solution = compressToMetaGraph(multi_graph, List_scaffolds)
    graph_to_graphml(meta_graph_solution, workdir, "solution_compressed")

    return meta_graph_solution


def delete_all_nodes_smaller_than(graph_ml, kmer_size, overlaps_min):
    workdir = os.path.dirname(os.path.abspath(str(graph_ml)))
    graph = nx.read_graphml(graph_ml)
    graph_copy = graph.copy()
    list_nodes_to_remove = []
    for node in graph.nodes():
        if graph.node[node]['UnitigLength'] < 2 * kmer_size - overlaps_min:
            list_nodes_to_remove.append(node)

    graph.remove_nodes_from(list_nodes_to_remove)
    file = os.path.basename(graph_ml)
    name = file.split('.')[0]
    graph_file = graph_to_graphml(graph, workdir, name)
    graph_to_graphml(graph_copy, workdir, name + "before_deletion")
    return graph_file


def delete_all_nodes_smaller_than_size(graph_ml, size_min):
    workdir = os.path.dirname(os.path.abspath(str(graph_ml)))
    graph = nx.read_graphml(graph_ml)
    graph_copy = graph.copy()
    list_nodes_to_remove = []
    for node in graph.nodes():
        if graph.node[node]['UnitigLength'] < size_min:
            list_nodes_to_remove.append(node)

    graph.remove_nodes_from(list_nodes_to_remove)
    file = os.path.basename(graph_ml)
    name = file.split('.')[0]
    graph_file = graph_to_graphml(graph, workdir, name)
    graph_to_graphml(graph_copy, workdir, name + "before_deletion")
    return graph_file


def remove_mp_between_smalls(graph):
    list_mp = [(u, v) for (u, v) in graph.edges() if graph.edge[u][v]['Distance'] >= 0]
    list_mp_to_remove = [(u, v) for (u, v) in list_mp if
                         graph.node[u]['UnitigLength'] < 250 or graph.node[v]['UnitigLength'] < 250]
    return list_mp_to_remove


def is_a_repeat(unitig_name):
    return unitig_name.split('_')[1] != '0'


def get_mp_graph_no_repeat(mp_graph):
    mp_graph_no_repeat = mp_graph.copy()
    links_repeat = [(u, v, d) for u, v, d in mp_graph.edges(data=True) if
                    d['Distance'] >= 0 and (is_a_repeat(u) or is_a_repeat(v))]
    mp_graph_no_repeat.remove_edges_from(links_repeat)
    return mp_graph_no_repeat


def get_mp_graph(graph_problem):
    mp_graph = graph_problem.copy()
    overlaps = [(u, v, d) for u, v, d in graph_problem.edges(data=True) if d['Distance'] < 0]
    mp_graph.remove_edges_from(overlaps)
    singletons = [node for node in mp_graph.nodes() if
                  len(mp_graph.successors(node)) == 0 and len(mp_graph.predecessors(node)) == 0]
    mp_graph.remove_nodes_from(singletons)
    return mp_graph


def reverse_graph_solution(graph_solution):
    reverse_graph = graph_solution.reverse()
    mapping = {}
    for node in graph_solution.nodes():
        mapping[node] = get_reverse(node)
    nx.relabel_nodes(reverse_graph, mapping)
    return reverse_graph


def map_to_graphml_with_mp(file_map, pb_file, workdir):
    dict_occ = {}
    list_nodes = []
    list_edges = []
    graph_pb = nx.read_graphml(pb_file)

    name = file_map.split('.')[0]

    with open(file_map) as f:
        lines = f.readlines()
        for line in lines:
            parts = line.split()
            orientation = 'R' if parts[2] == '-' else 'F'
            name_initial_unitig = parts[3]
            parts_name_initial = name_initial_unitig.split('__')
            number_unitig = parts_name_initial[0]
            if number_unitig in dict_occ:
                dict_occ[number_unitig] += 1
            else:
                dict_occ[number_unitig] = 0
            list_nodes.append(number_unitig + '_' + str(dict_occ[number_unitig]) + '_' + orientation)
    if len(list_nodes) >= 2 and list_nodes[0].split('_')[0] == list_nodes[-1].split('_')[0]:
        list_nodes.pop()
    for first, second in zip(list_nodes, list_nodes[1:]):
        list_edges.append((first, second))

    G = graph_pb.copy()
    print('list Edges : ', list_edges)
    for edge in list_edges:
        if edge in G.edges():
            print("OK")
        else:
            print('ERROR', edge)
    nodes_to_remove = [node for node in graph_pb.nodes() if node not in list_nodes]
    G.remove_nodes_from(nodes_to_remove)
    overlaps = [(source, sink) for (source, sink) in graph_pb.edges() if graph_pb.edge[source][sink]['Distance'] < 0]
    edges_to_remove = [(source, sink) for (source, sink) in overlaps if (source, sink) not in list_edges]

    G.remove_edges_from(edges_to_remove)

    list_good_links = [(source, sink) for (source, sink) in G.edges() if
                       G.edge[source][sink]['Distance'] >= 0 and (
                           G.node[source]['BigUnitig'] or G.node[sink]['BigUnitig'])]
    bad_links = [(source, sink) for (source, sink) in G.edges() if
                 G.edge[source][sink]['Distance'] >= 0 and (source, sink) not in list_good_links]

    G.remove_edges_from(remove_mp_between_smalls(G))
    G.remove_edges_from(bad_links)
    G.remove_nodes_from(
        [node for node in G.nodes() if not G.node[node]['UnitigLength'] >= 156])

    graph_to_graphml(G, workdir, name + '_solution')
    return G


def map_to_graphml(file_map, workdir):
    dict_occ = {}
    list_nodes = []
    list_edges = []

    name = file_map.split('.')[0]

    with open(file_map) as f:
        lines = f.readlines()
        for line in lines:
            parts = line.split()
            orientation = 'F' if parts[2] == '+' else 'R'
            name_initial_unitig = parts[3]
            parts_name_initial = name_initial_unitig.split('__')
            number_unitig = parts_name_initial[0]
            if number_unitig in dict_occ:
                dict_occ[number_unitig] += 1
            else:
                dict_occ[number_unitig] = 0
            list_nodes.append(number_unitig + '_' + str(dict_occ[number_unitig]) + '_' + orientation)
    if len(list_nodes) >= 2 and list_nodes[0].split('_')[0] == list_nodes[-1].split('_')[0]:
        list_nodes.pop()
    for first, second in zip(list_nodes, list_nodes[1:]):
        list_edges.append((first, second))

    G = nx.DiGraph()
    G.add_edges_from(list_edges)

    graph_to_graphml(G, workdir, name + '_solution')
    return G, list_nodes, name


def delete_singletons(graph):
    graph_copy = graph.copy()
    graph_copy.remove_nodes_from(nx.isolates(graph))
    return graph_copy


def create_graph_problem_from_step(file_step, name_subject, workdir):
    workdir = os.path.dirname(os.path.abspath(str(file_step)))

    workdir += "/Graphes"
    if not os.path.isdir(workdir):
        os.system("mkdir " + workdir)
    workdir += "/Problem"
    if not os.path.isdir(workdir):
        os.system("mkdir " + workdir)

    print(workdir)

    parsed_data = parser.parse_file(file_step)
    greatest_edge_length = max([edge.distance_max() for edge in parsed_data[1]])
    initGraph = create_graph_for_model(parsed_data[0], parsed_data[1], greatest_edge_length)
    # initGraph = delete_singletons(initGraph)
    assert len(initGraph.nodes()) != 0, "Erreur : create_graph graphe vide"
    graph_file = graph_to_graphml(initGraph, workdir, name_subject + "_init")
    return graph_file


def create_graphs_problem_from_step(file_step, name_subject):
    workdir = os.path.dirname(os.path.abspath(str(file_step)))

    workdir += "/Graphes"
    if not os.path.isdir(workdir):
        os.system("mkdir " + workdir)
    workdir += "/Problem"
    if not os.path.isdir(workdir):
        os.system("mkdir " + workdir)

    Nodes, Overlaps, Links = parser.parse_file_with_edges_division(file_step)

    if not Links:
        longest_link = 0
    else:
        longest_link = max([edge.distance_max() for edge in Links])
    OverlapsGraph = create_graph_for_model(Nodes, Overlaps, longest_link)
    #OverlapsGraph = clean_graph_from_source_and_sink(OverlapsGraph)
    OverlapsGraph = delete_singletons(OverlapsGraph)
    graph_file_overlaps = graph_to_graphml(OverlapsGraph, workdir, name_subject + '_overlaps')
    LinksGraph = create_graph_for_model(Nodes, Links, longest_link)
    # LinksGraph = clean_graph_from_source_and_sink(LinksGraph)
    # LinksGraph = delete_singletons(LinksGraph)
    graph_file_links = graph_to_graphml(LinksGraph, workdir, name_subject + '_links')

    return graph_file_overlaps, graph_file_links


def get_reverse(node):
    node_split = node.split('_')
    orientation = node_split[2]
    if orientation == "R":
        reverse = '' + node_split[0] + '_' + node_split[1] + '_' + 'F'
    else:
        reverse = '' + node_split[0] + '_' + node_split[1] + '_' + 'R'
    return reverse


def get_reverse_bis(node_without_occ):
    node_split = node_without_occ.split('__')
    orientation = node_split[1]
    if orientation == "R":
        reverse = '' + node_split[0] + '__' + 'F'
    else:
        reverse = '' + node_split[0] + '__' + 'R'
    return reverse


def delete_useless_sources_sinks(graph, list_source_sink_ok):
    list_source_sink_not_ok = [node for node in graph.nodes() if
                               node not in list_source_sink_ok and (
                                   graph.successors(node) == [] or graph.predecessors(node) == [])]

    while list_source_sink_not_ok:
        graph.remove_nodes_from(list_source_sink_not_ok)
        list_source_sink_not_ok = [node for node in graph.nodes() if
                                   node not in list_source_sink_ok and (
                                       graph.successors(node) == [] or graph.predecessors(node) == [])]


def multiply_nodes(nodes):
    dico_nodes_oriented = {}
    dico_nodes_length = {}
    list_nodes_for_graph = []
    for node in nodes:
        node_index_str = str(node.get_index())
        dico_nodes_oriented[node_index_str + "F"] = []
        dico_nodes_oriented[node_index_str + "R"] = []
        for i in range(node.get_count()):
            dico_nodes_oriented[node_index_str + "F"].append(str(node.get_index()) + "_" + str(i) + "_" + "F")
            list_nodes_for_graph.append(node_index_str + "_" + str(i) + "_" + "F")
            dico_nodes_length[node_index_str + "_" + str(i) + "_" + "F"] = node.get_size_total()
            dico_nodes_oriented[node_index_str + "R"].append(str(node.get_index()) + "_" + str(i) + "_" + "R")
            list_nodes_for_graph.append(node_index_str + "_" + str(i) + "_" + "R")
            dico_nodes_length[node_index_str + "_" + str(i) + "_" + "R"] = node.get_size_total()

    return dico_nodes_oriented, list_nodes_for_graph, dico_nodes_length


def make_edge_for_graph(dico_nodes_oriented, edges):
    list_edges_for_graph = []
    dico_edge_length = {}
    for edge in edges:
        source_str = str(edge.get_source()) + graph_tools.ORIENTATIONS_STR[edge.orientation_source()]
        sink_str = str(edge.get_sink()) + graph_tools.ORIENTATIONS_STR[edge.orientation_sink()]
        weight = edge.distance()
        weight_min = edge.distance_min()
        weight_max = edge.distance_max()
        for (source, sink) in itertools.product(dico_nodes_oriented[source_str], dico_nodes_oriented[sink_str]):
            if source != sink:
                list_edges_for_graph.append((source, sink))
                dico_edge_length[(source, sink)] = (weight, weight_min, weight_max)
    return list_edges_for_graph, dico_edge_length


def create_graph_for_model(nodes, edges, lenght_min_of_a_big_unitig):
    G = nx.DiGraph()
    dico_nodes_oriented, list_nodes_for_graph, dico_nodes_length = multiply_nodes(nodes)
    list_edges_for_graph, dico_edge_length = make_edge_for_graph(dico_nodes_oriented, edges)
    for node in list_nodes_for_graph:
        if dico_nodes_length[node] > lenght_min_of_a_big_unitig:
            G.add_node(node, UnitigLength=dico_nodes_length[node], BigUnitig=True, Fusion=False)
        else:
            G.add_node(node, UnitigLength=dico_nodes_length[node], BigUnitig=False, Fusion=False)
    for (edge_source, edge_sink) in list_edges_for_graph:
        if dico_edge_length[(edge_source, edge_sink)][0] > 0:
            if G.node[edge_source]["UnitigLength"] > LIMIT_SIZE and G.node[edge_sink]["UnitigLength"] > LIMIT_SIZE:
                G.add_edge(edge_source, edge_sink, Distance=dico_edge_length[(edge_source, edge_sink)][0],
                           Distance_min=dico_edge_length[(edge_source, edge_sink)][1],
                           Distance_max=dico_edge_length[(edge_source, edge_sink)][2])
        else:
            G.add_edge(edge_source, edge_sink, Distance=dico_edge_length[(edge_source, edge_sink)][0],
                       Distance_min=dico_edge_length[(edge_source, edge_sink)][1],
                       Distance_max=dico_edge_length[(edge_source, edge_sink)][2])
    return G


def graph_to_graphml(G, workdir, name):
    nx.draw(G)
    print(workdir + "/" + name + ".graphml")
    nx.write_graphml(G, workdir + "/" + name + ".graphml")
    return workdir + "/" + name + ".graphml"


def map_to_graph(map_name, workdir, name, graph_problem_initial):
    Graphe_forward = graph_problem_initial.copy()
    Graphe_reverse = graph_problem_initial.copy()

    file_map = open(map_name, "r")
    line = file_map.readline()
    List_line = []
    List_number_unitig = []

    while line != "":
        split = line.split()
        orientation = graph_tools.ORIENTATIONS_STR_to_STR[split[2]]
        List_line.append(split[3].split('__')[0] + '__' + orientation)
        List_number_unitig.append(split[3].split('__')[0])
        line = file_map.readline()

    unitigs_begin_end = List_line[0]
    unitigs_without_begin_end = List_line[1:-1]
    unitigs_occurence = {k: List_number_unitig.count(k) for k in set(List_number_unitig)}

    list_number_unitig_unique = [unitig_number for unitig_number in unitigs_occurence if
                                 unitigs_occurence[unitig_number] == 1]
    list_unitig_solution_forward = []

    for unitig in unitigs_without_begin_end:
        unitig_split = unitig.split('__')
        if unitig.split('__')[0] == unitigs_begin_end.split('__')[0]:
            list_unitig_solution_forward.append(
                unitig_split[0] + '_' + str(unitigs_occurence[unitig.split('__')[0]] - 2) + '_' + unitig_split[1])
            unitigs_occurence[unitig.split('__')[0]] -= 2
        else:
            list_unitig_solution_forward.append(
                unitig_split[0] + '_' + str(unitigs_occurence[unitig.split('__')[0]] - 1) + '_' + unitig_split[1])
            unitigs_occurence[unitig.split('__')[0]] -= 1

    unitigs_begin_end_split = unitigs_begin_end.split('__')
    unitigs_begin_end = unitigs_begin_end_split[0] + '_' + str(
        unitigs_occurence[unitigs_begin_end.split('__')[0]] - 1) + '_' + \
                        unitigs_begin_end_split[1]
    list_unitig_solution_forward = [unitigs_begin_end] + list_unitig_solution_forward + [unitigs_begin_end]
    list_unitig_solution_reverse = [get_reverse(node) for node in list_unitig_solution_forward]
    list_unitig_solution_reverse.reverse()


    list_edges_solution_forward = []
    iter_list = iter(list_unitig_solution_forward)
    source = next(iter_list)
    for i in range(len(list_unitig_solution_forward) - 1):
        sink = next(iter_list)
        list_edges_solution_forward.append((source, sink))
        source = sink

    list_edges_solution_reverse = []
    iter_list = iter(list_unitig_solution_reverse)
    source = next(iter_list)
    for i in range(len(list_unitig_solution_reverse) - 1):
        sink = next(iter_list)
        list_edges_solution_reverse.append((source, sink))
        source = sink

    for (source, sink) in graph_problem_initial.edges():
        if (source, sink) not in list_edges_solution_forward:
            Graphe_forward.remove_edge(source, sink)
        if (source, sink) not in list_edges_solution_reverse:
            Graphe_reverse.remove_edge(source, sink)

    for node in graph_problem_initial.nodes():
        if node not in list_unitig_solution_forward:
            Graphe_forward.remove_node(node)
        if node not in list_unitig_solution_reverse:
            Graphe_reverse.remove_node(node)
    print(unitigs_occurence)
    first_unique_unitig = [node for node in list_unitig_solution_forward if
                           node.split('_')[0] in list_number_unitig_unique]
    first_unique_unitig = [node for node in first_unique_unitig if graph_problem_initial.node[node]["BigUnitig"]][0]

    graph_to_graphml(Graphe_forward, workdir, name + "_forward")
    graph_to_graphml(Graphe_reverse, workdir, name + "_reverse")

    return Graphe_forward, Graphe_reverse, first_unique_unitig


def complete_graph_with(nodes, edges, graphe_solution, graphe_global):
    for node in nodes:
        if node not in graphe_solution.nodes():
            UnitLength = graphe_global.node[node]['UnitigLength']
            graphe_solution.add_node(node, UnitigLength=UnitLength, BigUnitig=graphe_global.node[node]['BigUnitig'],
                                     Fusion=False)

    for (source, sink) in edges:
        if (source, sink) in graphe_global.edges():
            d = graphe_global.edge[source][sink]['Distance']
            d_min = graphe_global.edge[source][sink]['Distance_min']
            d_max = graphe_global.edge[source][sink]['Distance_max']
            graphe_solution.add_edge(source, sink, Distance=d, Distance_min=d_min, Distance_max=d_max)


def get_edges_used_from_solution(solution):
    edges_used = []
    for element in solution:
        if "e_" in element:
            element_split = element.split('_')
            source = element_split[1] + '_' + element_split[2] + '_' + element_split[3]
            sink = element_split[4] + '_' + element_split[5] + '_' + element_split[6]
            edges_used.append((source, sink))
    return edges_used


def get_number_cycle(solution):
    G = nx.DiGraph()
    edges_used = get_edges_used_from_solution(solution)
    G.add_edges_from(edges_used)
    print("il y a " + str(len(list(nx.simple_cycles(G)))) + " cycles dans la solution")
    return len(list(nx.simple_cycles(G)))


def get_cycles(solution):
    G = nx.DiGraph()
    edges_used = get_edges_used_from_solution(solution)
    G.add_edges_from(edges_used)
    return list(nx.simple_cycles(G))


def split_solution_step1(solution):
    nodes = []
    edges = []
    for element in solution:
        if 'y_' in element or 'o_' in element or 'i_' in element:
            node = element[2:]
            nodes.append(node)
        if 'e_' in element:
            edge = element[2:]
            edge_split = edge.split('_')
            source = edge_split[0] + '_' + edge_split[1] + '_' + edge_split[2]
            sink = edge_split[3] + '_' + edge_split[4] + '_' + edge_split[5]
            edges.append((source, sink))
    return nodes, edges


def save_solution_step1_unique(list_solutions, graphe_problem, workdir, name_subject, number_solution):
    solution = list_solutions[0]
    G = nx.DiGraph()
    nodes, edges = split_solution_step1(solution)
    for node in nodes:
        if 'Fusion' in graphe_problem.node[node]:
            G.add_node(str(node), UnitigLength=graphe_problem.node[node]['UnitigLength'],
                       BigUnitig=graphe_problem.node[node]['BigUnitig'], Fusion=graphe_problem.node[node]['Fusion'])
        else:
            G.add_node(str(node), UnitigLength=graphe_problem.node[node]['UnitigLength'],
                       BigUnitig=graphe_problem.node[node]['BigUnitig'], Fusion=False)

    for (source, sink) in edges:
        G.add_edge(source, sink, Distance=graphe_problem.edge[source][sink]['Distance'],
                   Distance_min=graphe_problem.edge[source][sink]['Distance_min'],
                   Distance_max=graphe_problem.edge[source][sink]['Distance_max'])

    graph_to_graphml(G, workdir, name_subject + "_step1_solution" + str(number_solution))
    return G


def save_solution_step1(list_solutions, graphe_problem, workdir, name_subject):
    number_solution = 0
    for solution in list_solutions:
        G = nx.DiGraph()
        nodes, edges = split_solution_step1(solution)
        for node in nodes:
            G.add_node(str(node), UnitigLength=graphe_problem.node[node]['UnitigLength'], BigUnitig=True, Fusion=False)
        for (source, sink) in edges:
            G.add_edge(source, sink, Distance=graphe_problem.edge[source][sink]['Distance'],
                       Distance_min=graphe_problem.edge[source][sink]['Distance_min'],
                       Distance_max=graphe_problem.edge[source][sink]['Distance_max'])
        graph_to_graphml(G, workdir, name_subject + "_step1_solution" + str(number_solution))
        number_solution += 1


def create_graph_step1_with_split_nodes(workdir, graph_problem, node_to_split, free_index_for_split):
    graphe_copy = graph_problem.copy()
    graphe_copy.remove_node(node_to_split)
    succs = graph_problem.successors(node_to_split)
    preds = graph_problem.predecessors(node_to_split)
    node_source = str(free_index_for_split) + "_0_F"
    node_sink = str(free_index_for_split + 1) + "_0_F"
    graphe_copy.add_node(node_source, UnitigLength=graph_problem.node[node_to_split]['UnitigLength'], BigUnitig=True,
                         Fusion=False)
    graphe_copy.add_node(node_sink, UnitigLength=graph_problem.node[node_to_split]['UnitigLength'], BigUnitig=True,
                         Fusion=False)

    for succ in succs:
        graphe_copy.add_edge(node_source, succ, Distance=graph_problem.edge[node_to_split][succ]['Distance'],
                             Distance_min=graph_problem.edge[node_to_split][succ]['Distance_min'],
                             Distance_max=graph_problem.edge[node_to_split][succ]['Distance_max'])
    for pred in preds:
        graphe_copy.add_edge(pred, node_sink, Distance=graph_problem.edge[pred][node_to_split]['Distance'],
                             Distance_min=graph_problem.edge[pred][node_to_split]['Distance_min'],
                             Distance_max=graph_problem.edge[pred][node_to_split]['Distance_max'])
    graph_to_graphml(graphe_copy, workdir, "graph_problem_step1_after_split")
    return graphe_copy, node_source, node_sink


def create_graph_global(graph_solution_pred, graph_problem, list_source_sink_holes):
    for node in graph_solution_pred:
        reverse_node = get_reverse(node)
        if node not in list_source_sink_holes and node in graph_problem:
            graph_problem.remove_node(node)
        if reverse_node in graph_problem.nodes():
            graph_problem.remove_node(reverse_node)


def get_weakly_connected_component(graph):
    test = list(nx.weakly_connected_components(graph))
    return test, len(test)


def get_path_from_source_to_sink(graph, source, sink):
    print("source : ", source)
    print("sink : ", sink)
    list_path = list(nx.all_simple_paths(graph, source, sink))
    print(list_path)
    assert len(list_path) <= 1, 'Erreur : plusieurs chemin dans une meme composante connexe lineaire'
    return list_path


def get_path_from_source_to_sink_2(graph, source, sink):
    list_path = list(nx.all_simple_paths(graph, source, sink))
    if not list_path:
        return [[source]]
    print(list_path)
    assert len(list_path) <= 1, 'Erreur : plusieurs chemin dans une meme composante connexe lineaire'
    return list_path


def get_source_and_sink(graphe, list_nodes):
    return [node for node in list_nodes if graphe.successors(node) == []], [node for node in list_nodes if
                                                                            graphe.predecessors(node) == []]


def clean_graph_from_source_and_sink(graphe):
    list_nodes = graphe.nodes()
    nodes_to_delete = [node for node in list_nodes if graphe.successors(node) == [] or graphe.predecessors(node) == []]
    while (len(nodes_to_delete) != 0):
        graphe.remove_nodes_from(nodes_to_delete)
        list_nodes = graphe.nodes()
        nodes_to_delete = [node for node in list_nodes if
                           graphe.successors(node) == [] or graphe.predecessors(node) == []]
    return graphe


def add_nodes_and_overlaps_to_multigraph(multigraph, ovl_graph):
    nodes = ovl_graph.nodes()
    for node in nodes:
        multigraph.add_node(node, UnitigLength=ovl_graph.node[node]['UnitigLength'],
                            BigUnitig=ovl_graph.node[node]['BigUnitig'])
    for (u, v) in ovl_graph.edges():
        multigraph.add_edge(u, v, Distance=ovl_graph.edge[u][v]['Distance'], Type="overlaps")


def add_links_to_multigraph(multigraph, links_graph):
    for (u, v) in links_graph.edges():
        multigraph.add_edge(u, v, Distance=links_graph.edge[u][v]['Distance'],
                            Distance_min=links_graph.edge[u][v]['Distance_min'],
                            Distance_max=links_graph.edge[u][v]['Distance_max'], Type="links")


def delete_all_singletons(multigraph):
    multigraph.remove_nodes_from(nx.isolates(multigraph))
    return multigraph


def add_links_to_multigraph_if_repetition_factor_less_than(multigraph, links_graph, repetition_factor):
    for (u, v) in links_graph.edges():
        u_split = u.split('_')
        v_split = v.split('_')
        if int(u_split[1]) <= repetition_factor and int(v_split[1]) <= repetition_factor:
            multigraph.add_edge(u, v, Distance=links_graph.edge[u][v]['Distance'],
                                Distance_min=links_graph.edge[u][v]['Distance_min'],
                                Distance_max=links_graph.edge[u][v]['Distance_max'], Type="links")


def create_multigraph_problem(ovl_graph_file, links_graph_file, workdir):
    ovl_graph = nx.read_graphml(ovl_graph_file)
    links_graph = nx.read_graphml(links_graph_file)

    multigraph = nx.MultiDiGraph()

    add_nodes_and_overlaps_to_multigraph(multigraph, ovl_graph)
    add_links_to_multigraph(multigraph, links_graph)
    name = os.path.basename(ovl_graph_file).split('_')[0]
    # multigraph = delete_all_singletons(multigraph)
    # gm.clean_graph_from_source_and_sink(multigraph)
    nx.write_graphml(multigraph, workdir+'/' + name + '_multigraph.graphml')


def create_multigraph_problem_simulator(ovl_graph_file, links_graph_file, workdir):
    name = os.path.basename(ovl_graph_file).split('_')[0]
    ovl_graph = nx.read_graphml(ovl_graph_file)
    links_graph = nx.read_graphml(links_graph_file)

    multigraph = nx.MultiDiGraph()
    nx.write_graphml(ovl_graph, workdir +'/' + name + '_debugOverlaps.graphml')
    add_nodes_and_overlaps_to_multigraph(multigraph, ovl_graph)
    add_links_to_multigraph(multigraph, links_graph)

    print("MultiGraph creation :")
    print(workdir + name + '_multigraph.graphml')
    nx.write_graphml(multigraph, workdir + name + '_multigraph.graphml')
    return multigraph

def diGraphToMultiDiGraph(mg_graph):
    if mg_graph.is_multigraph():
        return mg_graph
    else:
        multigraph_mg_graph = nx.MultiDiGraph()
        for node in mg_graph.nodes():
            multigraph_mg_graph.add_node(node, UnitigLength=mg_graph.node[node]['UnitigLength'],
                                         BigUnitig=mg_graph.node[node]['BigUnitig'])
        for (u, v) in mg_graph.edges():
            multigraph_mg_graph.add_edge(u, v, Distance=mg_graph.edge[u][v]['Distance'], Type="overlaps")
        return multigraph_mg_graph


def build_DCEP_multigraph(mg_graph):
    mg_graph = diGraphToMultiDiGraph(mg_graph)
    multigraph_DCEP = nx.MultiDiGraph()
    for node in mg_graph.nodes():
        multigraph_DCEP.add_node(node, BigUnitig=mg_graph.node[node]['BigUnitig'])
    for (u, v, k) in mg_graph.edges(keys=True):
        if mg_graph.edge[u][v][k]['Type'] == 'overlaps':
            multigraph_DCEP.add_edge(u, v,
                                     Distance=mg_graph.edge[u][v][k]['Distance'] + mg_graph.node[u]['UnitigLength'],
                                     Type="overlaps")
        elif mg_graph.edge[u][v][k]['Type'] == 'links':
            multigraph_DCEP.add_edge(u, v,
                                     Distance=mg_graph.edge[u][v][k]['Distance'] + mg_graph.node[u]['UnitigLength'],
                                     Distance_min=mg_graph.edge[u][v][k]['Distance_min'] + mg_graph.node[u][
                                         'UnitigLength'],
                                     Distance_max=mg_graph.edge[u][v][k]['Distance_max'] + mg_graph.node[u][
                                         'UnitigLength'], Type="links")
    return multigraph_DCEP


def make_dico_node(mg_graph):
    dico_nodes = {}
    for node in mg_graph.nodes():
        dico_nodes[node] = mg_graph.node[node]['UnitigLength']
    return dico_nodes


if __name__ == '__main__':
    G = nx.DiGraph()
    print(G.nodes())
    print(G.edges())
    nx.draw(G)
    nx.write_graphml(G, "test.graphml")
