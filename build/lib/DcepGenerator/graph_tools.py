#!/usr/bin/env python

# -*- coding: utf-8 -*-
"""
Created on Fri May 29 10:02:57 2015

@author: seb
"""
ORIENTATIONS_STR_to_STR = {'+': 'F', '-': 'R'}
ORIENTATIONS_STR = {True: 'F', False: 'R'}
INV_ORIENTATION = {'F': True, '+': True, 'R': False, '-': False}


class Solution:
    def __init__(self, coord, genome, contig_index, size, coord2, contig,
                 orientation, overlapping):
        self._coord = coord
        self._genome = genome
        self._contig_index = contig_index
        self._size = size
        self._coord2 = coord2
        self._contig = contig
        self._orientation = orientation
        self._overlapping = overlapping

    def __repr__(self):
        return "<Solution: coord={}, genome={}, contig_index={}, size={}, coord2={}, contig={}, orientation={}, overlapping={}>".format(
            self._coord, self._genome, self._contig_index, self._size, self._coord2, self._contig, self._orientation,
            self._overlapping)

    def shortRepr(self, nodes):
        return "{:<15} {}".format(nodes[self._contig_index].name(),
                                  ORIENTATIONS_STR[self._orientation])


class Edge:
    def __init__(self, source, sink, dist, orient_source, orient_sink, distmin=0, distmax=0):
        self._source = source
        self._sink = sink
        self._dist = dist
        if distmin == 0:
            self._distmin = dist
        else:
            self._distmin = distmin
        if distmax == 0:
            self._distmax = dist
        else:
            self._distmax = distmax
        self._orient_source = orient_source
        self._orient_sink = orient_sink

    def __hash__(self):
        return hash((self._source, self._sink, self._dist, self._orient_source, self._orient_sink))

    def __eq__(self, other):
        return (self._dist, self._orient_sink, self._orient_source) == (
            other._dist, other._orient_sink, other._orient_source)

    def __repr__(self):
        return "<Edge: C{}{}, C{}{}, dist={}, dist_min={}, dist_max={}>".format(
            self._source, ORIENTATIONS_STR[self._orient_source], self._sink, ORIENTATIONS_STR[self._orient_sink],
            self._dist, self._distmin, self._distmax)

    def shortRepr(self, nodes, formating_str="{:<15} {:<15} {} {} {:>6}",
                  format_funciton=None):
        def fields_reformatter(source, sink, orient_source, orient_sink, dist):
            return (nodes[source].name(), nodes[sink].name(),
                    ORIENTATIONS_STR[orient_source], ORIENTATIONS_STR[orient_sink], dist)

        format_funciton = format_funciton or fields_reformatter

        return formating_str.format(*format_funciton(self._source, self._sink,
                                                     self._orient_source, self._orient_sink, self._dist))

    def distance(self):
        return self._dist

    def distance_max(self):
        return self._distmax

    def distance_min(self):
        return self._distmin

    def set_distance(self, new_distance):
        self._dist = new_distance

    def reverse(self):
        return Edge(self._sink, self._source, self._dist, not self._orient_sink,
                    not self._orient_source)

    def get_symetrique(self):
        return Edge(self._sink, self._source, self._dist, not self._orient_sink, not self._orient_source)

    def get_source(self):
        return self._source

    def get_sink(self):
        return self._sink

    def orientation_source(self):
        return self._orient_source

    def orientation_sink(self):
        return self._orient_sink

    def get_contig_oriented_source(self):
        return Contig_Oriented(self._source, self._orient_source)

    def get_contig_oriented_sink(self):
        return Contig_Oriented(self._sink, self._orient_sink)

    def contain_Node(self, node):
        #      print(str(self._source))
        #      print(str(node.get_index()))
        #      print(str(self._sink))
        #      print(str(node.get_index()))
        test = self._source == node.get_index() or self._sink == node.get_index()
        #      print(test)
        return test

    def getOtherOrientation(self, node_number, orientation):
        if node_number == self._source:
            if self._orient_source == orientation:
                return self._orient_sink
            else:
                return not self._orient_sink
        else:
            if self._orient_sink == orientation:
                return self._orient_source
            else:
                return not self._orient_source


class Contig_Oriented:
    def __init__(self, index, orientation):
        self._index = index
        self._orientation = orientation

    def __hash__(self):
        return hash((self._index, self._orientation))

    def __eq__(self, other):
        return (self._index, self._orientation) == (
            other._index, other._orientation)

    def __repr__(self):
        return "<Contig: C{}{}>".format(
            self._index, ORIENTATIONS_STR[self._orientation])

    def get_index(self):
        return self._index

    def get_orientation(self):
        return self._orientation

    def get_reverse(self):
        return Contig_Oriented(self._index, not self._orientation)


class Node:
    def __init__(self, index, size, count):
        self._index = index
        self._size = size
        self._size_right = size // 2
        self._size_left = self._size - self._size_right
        self._count = count

    def __hash__(self):
        return hash((self._index, self._size, self._count))

    def __eq__(self, other):
        return (self._index, self._size, self._count) == (
            other._index, other._size, other._count)

    def __repr__(self):
        return "<Node: C{}, size={}, count={}>".format(
            self._index, self._size, self._count)

    def give_me_the_forward_version(self):
        return Contig_Oriented(self._index, True)

    def give_me_the_reverse_version(self):
        return Contig_Oriented(self._index, False)

    def get_count(self):
        return self._count

    def decrease_count_by_one(self):
        self._count = self._count - 1

    def increase_count_by_one(self):
        self._count = self._count + 1

    def get_index(self):
        return self._index

    def get_size_total(self):
        return self._size

    def get_size_right(self):
        return self._size_right

    def get_size_left(self):
        return self._size_left

    def get_clone(self):
        return Node(self._index, self._size, 0)

    # for index, node in enumerate(nodes)
    def similar(self, other):
        return (self._size, self._count) == (other._size, other._count)

    def name(self):
        return "{}__len__{}".format(self._index, self._size)

    def shortRepr(self, _, formating_str="{:<15} {:<8} {}", format_funciton=None):
        def fields_reformatter(index, size, count):
            return (self.name(), size, count)

        format_funciton = format_funciton or fields_reformatter

        return formating_str.format(*format_funciton(self._index, self._size,
                                                     self._count))


class Link:
    def __init__(self, contig_source, contig_orientation_source, contig_sink, contig_orientation_sink, distance_min,
                 distance_max):
        self._contig_source = contig_source
        self._contig_orientation_source = contig_orientation_source
        self._contig_sink = contig_sink
        self._contig_orientation_sink = contig_orientation_sink
        self._distance_min = distance_min
        self._distance_max = distance_max
        self._distance_average = (self._distance_min + self._distance_max) // 2

    def __hash__(self):
        return hash((self._contig_source, self._contig_orientation_source,
                     self._contig_sink, self._contig_orientation_sink,
                     self._distance_min, self._distance_max))

    def __eq__(self, other):
        return ((self._contig_source, self._contig_orientation_source,
                 self._contig_sink, self._contig_orientation_sink,
                 self._distance_min, self._distance_max)
                ==
                (other._contig_source, other._contig_orientation_source,
                 other._contig_sink, other._contig_orientation_sink,
                 other._distance_min, other._distance_max))

    def __str__(self):
        return "<Link: contig_source= {}, contig_orientation_source= {},contig_sink= {}, contig_orientation_sink= {}, distance_min= {}, distance_max= {}, distance_average= {}>".format(
            self._contig_source, self._contig_orientation_source, self._contig_sink, self._contig_orientation_sink,
            self._distance_min, self._distance_max, self._distance_average)

    def lower_bound(self):
        return self._distance_min

    def upper_bound(self):
        return self._distance_min

    #    a corriger !
    def distance_average(self):
        return self._distance_min

    def get_symetrique(self):
        return Link(self._contig_sink, not self._contig_orientation_sink, self._contig_source,
                    not self._contig_orientation_source, self._distance_min, self._distance_max)

    def get_edge(self):
        return Edge(self._contig_source, self._contig_sink, self._distance_min, self._contig_orientation_source,
                    self._contig_orientation_sink, self._distance_min, self._distance_max)
