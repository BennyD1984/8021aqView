
# !/usr/bin/env python
# -*- coding: utf-8 -*-
# ===================================================
# 8021aqView_v2.0 - created by Benjamin Domroese in 2021.
#
# Improvements planned: - SPSource ID Name display in MCast.
#                       - ISID view (number of connections per edge,
#                         not only by BVLAN usage, also by real ISID usage.
#
# Known Issues: Split SPB Domain not supported!

# CHANGELOG
# v2_0: Added APSP ECT counter per BVLAN.
# v2_0: Added "max_bvl_to_test" and "low_bvl_to_test", to reduce simulation
#       time in large networks. Change default range to 4000 - 4015.
# v2_0: Added "d_node_xor", to precalculate xor-values for all nodes
#       and all bvlans.
# v2_0: Added "d_tiebreaker_ect", for statistics on which byte a tiebreaker
#       was found for ECT.
# v2_0: Changed ECT algorithm to mechanism mentioned in:
#       "aq-ashwood-ECT-framework-1109-v1.pptx".
# v2_0: Added "d_edges_weight", for displaying loaded edges with metric.
# v2_0: Change ECT path counter to not include zero hop path (node to itself)

import math
import time
import os
from itertools import product
from prettytable import PrettyTable
from spbnetworks import d_node


class AbstractDijkstraSPF():

    def __init__(self, G, s, bvl):
        """ Calculate shortest path from s to other nodes in G. """
        self.__dist = dist = dict()
        self.__prev = prev = dict()
        self.ectcnt = 0
        visited = set()
        queue = set()

        dist[s] = 0
        prev[s] = s
        queue.add(s)

        while queue:
            u = min(queue, key=dist.get)
            for v in self.get_adjacent_nodes(G, u):
                if v in visited:
                    continue
                alt = self.get_distance(u) + self.get_edge_weight(G, u, v)

                # Check if alternate path is cheaper than actual path.
                if alt < self.get_distance(v):
                    dist[v] = alt
                    prev[v] = u
                    queue.add(v)

                # Check for equal cost path AND less hops on alternate path.
                elif alt == self.get_distance(v) and \
                        len(self.get_path(v))-1 > len(self.get_path(u)):
                    dist[v] = alt
                    prev[v] = u
                    queue.add(v)

                # Check for equal cost path AND equal number of hops.
                elif alt == self.get_distance(v) and \
                        len(self.get_path(v))-1 == len(self.get_path(u)):

                    # Get hop (used) before v, that be taken for XOR
                    # calculation.
                    v2 = self.get_path(v)[-2]

                    # v2_0: Get PATHID and find lowest SystemID
                    # (after XOR calc).
                    # Get all unique Nodes between FORK and JOIN
                    # points of subpath_used and subpath_alt.
                    subpath_used = list(filter(lambda x: x not in self.
                                        get_path(u), self.get_path(v2)))
                    subpath_alt = list(filter(lambda x: x not in self.
                                              get_path(v2), self.get_path(u)))

                    # Get LowPathID for subpath_used/alt.
                    # This is done by using the lowest first byte of (xored)
                    # SystemIDs in subpath_used, and compare it with the
                    # lowest first byte of (xored) SystemIDs in subpath_alt.
                    for i in range(8):
                        list_bytes_subpath_used = []
                        list_bytes_subpath_alt = []
                        for hop in subpath_used:
                            list_bytes_subpath_used.append(
                                [d_node_xor[int(bvl)][hop][i]])
                            xor_v = min(list_bytes_subpath_used)
                        for hop in subpath_alt:
                            list_bytes_subpath_alt.append(
                                [d_node_xor[int(bvl)][hop][i]])
                            xor_u = min(list_bytes_subpath_alt)
                        self.ectcnt += 1
                    # If min of first byte is the same, continue until byte 8.
                    # If no difference would be found, a duplicate SystemID
                    # would be present which causes an exception.
                        if xor_v != xor_u:
                            # For statistic purpose, save tiebreaker position.
                            d_tiebreaker_byte[i] = d_tiebreaker_byte[i]+1
                            break

                    # Compare XORed values.
                    if xor_v < xor_u:
                        dist[v] = alt
                        prev[v] = v2
                        queue.add(v)
                    elif xor_v > xor_u:
                        prev[v] = u

            queue.remove(u)
            visited.add(u)

    @staticmethod
    def get_adjacent_nodes(G, u):
        raise NotImplementedError()

    @staticmethod
    def get_edge_weight(G, u, v):
        raise NotImplementedError()

    def get_distance(self, u):
        """ Return the length of shortest path from s to u. """
        return self.__dist.get(u, math.inf)

    def get_path(self, v):
        """ Return the shortest path to v. """
        path = [v]
        while self.__prev[v] != v:
            v = self.__prev[v]
            path.append(v)
        return path[::-1]


class DijkstraSPF(AbstractDijkstraSPF):

    @staticmethod
    def get_adjacent_nodes(G, u):
        return G.get_adjacent_nodes(u)

    @staticmethod
    def get_edge_weight(G, u, v):
        return G.get_edge_weight(u, v)

    @staticmethod
    def byte_xor(self, ba1, ba2):
        return bytes([a ^ b for a, b in zip(ba1, ba2)])


class Graph():

    def __init__(self, adjacency_list=dict(), edge_weights=dict()):
        self.__adjacency_list = adjacency_list.copy()
        self.__edge_weights = edge_weights.copy()

    def add_edge(self, u, v, w):
        """ Add a new edge u -> v to graph with edge weight w. """
        self.__edge_weights[u, v] = w
        if u not in self.__adjacency_list:
            self.__adjacency_list[u] = set()
        self.__adjacency_list[u].add(v)

    def get_edge_weight(self, u, v):
        """ Get edge weight of edge between u and v. """
        return self.__edge_weights[u, v]

    def get_adjacent_nodes(self, u):
        """ Get nodes adjacent to u. """
        return self.__adjacency_list.get(u, set())

    def get_number_of_nodes(self):
        """ Return the total number of nodes in graph. """
        return len(self.__adjacency_list)

    def get_nodes(self):
        """ Return all nodes in this graph. """
        return self.__adjacency_list.keys()


def byte_xor(ba1, ba2):
    return bytes([a ^ b for a, b in zip(ba1, ba2)])


def get_max_possible_path_a_to_b(a, b):
    set_unique_path = set()  # Remove duplicate paths, keep only unique.
    for bvl in range(low_bvl_to_test, max_bvl_to_test):
        path_as_string = str(d_dijkstra[a][bvl].get_path(b))
        set_unique_path.add(path_as_string)  # Add string to set
    return (set_unique_path)


def get_max_possible_path_a_to_b_filtered(a, b, listbvls):
    set_unique_path = set()  # Remove duplicate paths, keep only unique.
    for bvl in range(low_bvl_to_test, max_bvl_to_test):
        if bvl in listbvls:
            path_as_string = str(d_dijkstra[a][bvl].get_path(b))
            set_unique_path.add(path_as_string)  # Add string to set.
    return (set_unique_path)


def get_max_possible_path_table():
    for k in sorted(d_ect_path, reverse=True):
        print("Number of connections with maximum {:>3d} possible ECT pathes: \
         {:>4d}".format(int(k), int(len(d_ect_path[k]))))
    return()


def get_path_a_to_b(a, b):
    print("Possible pathes from {} to {}:".format(a, b))
    print("===========================================")
    for bvl in range(low_bvl_to_test, max_bvl_to_test):
        print("BVLAN {}: {}".format(bvl,
              " -> ".join(d_dijkstra[a][bvl].get_path(b))))


def get_path_a_to_b_filtered(a, b, listbvls):
    print("Pathes from {} to {} in used BVLAN`s:".format(a, b))
    print("===========================================")
    for bvl in range(low_bvl_to_test, max_bvl_to_test):
        if bvl in listbvls:
            print("BVLAN {}: {}".format(bvl,
                  " -> ".join(d_dijkstra[a][bvl].get_path(b))))


def get_path_a_to_b_mc(a, b, bvl):
    print("BVLAN {}: {}".format(bvl, d_dijkstra[a][bvl].get_path(b)))


def main_menu():
    print(r" ___   __  ___     _                 __   __ _")
    print(r"( _ ) /  \|_  )   / | __ _  __ _  ___\ \ / /(_) ___ __ __ __")
    print(r"/ _ \| () |/ /  _ | |/ _` |/ _` ||___|\ V / | |/ -_)\ V  V /")
    print(r"\___/ \__//___|(_)|_|\__,_|\__, |      \_/  |_|\___| \_/\_/")
    print(r"                              |_|")
    print(" 1: Show Node List")
    print(" 2: Maximum possible path diversity as summary")
    print(" 3: Maximum possible path diversity between A and B")
    print(" 4: Display connections with n possible pathes")
    print("   ---------- Multicast view ---------------")
    print("10: Number of Multicast states per Node")
    print("11: Multicast states of a specific node")
    print("   ---------- Edge view --------------------")
    print("20: Show all edges")
    print("21: Usage of edges as summary")
    print("22: Usage of a specific edge")
    print("   ---------- Performance view -------------")
    print("30: Number of ECT Calculations")
    print("31: Time for calculation of SPF and MCast")
    print("   -----------------------------------------")
    print(" 0: Quit program\n")

    # Interactive usermenu.
    userchoice = input("Select Option:")
    if userchoice == "1":
        menu1()
    elif userchoice == "2":
        menu2()
    elif userchoice == "3":
        menu3()
    elif userchoice == "4":
        menu4()
    elif userchoice == "10":
        menu10()
    elif userchoice == "11":
        menu11()
    elif userchoice == "20":
        menu20()
    elif userchoice == "21":
        menu21()
    elif userchoice == "22":
        menu22()
    elif userchoice == "30":
        menu30()
    elif userchoice == "31":
        menu31()
    elif userchoice == "0":
        menu_exit()
    else:
        print("Wrong input! Please choose valid option.")
        main_menu()


def key_error_node():
    input1 = input("Key Error! Please check node name. " +
                   "Display list of nodes (y/n)?:")
    if input1 == "y" or input1 == "Y":
        print("List of Nodes: ")
        print("============== ")
        for k in d_node.keys():
            print(k)


def decorator_cls_and_return_mainmenu(func):
    def wrapper():
        os.system("cls")
        func()
        input("\nType anykey to continue")
        os.system("cls")
        main_menu()
    return wrapper


@decorator_cls_and_return_mainmenu
def menu1():
    print("List of nodes: ")
    print("============== ")
    for k in d_node.keys():
        print(k)


@decorator_cls_and_return_mainmenu
def menu2():
    print("Connections are counted bidirectional (A -> B and B -> A)")
    print("Counted are all paths available in BVLAN range {} to {}".format(
        low_bvl_to_test, max_bvl_to_test - 1))
    print("=========================================================")
    get_max_possible_path_table()


@decorator_cls_and_return_mainmenu
def menu3():
    input1 = input("Select source node: ")
    input2 = input("Select destination node: ")
    input3 = input("Only for used BVLAN`s (y/n): ")
    print("===========================================")
    try:
        if input3 == "y" or input3 == "Y":
            get_path_a_to_b_filtered(input1, input2, set_bvls_in_use)
            print("\nNumber of different pathes available " +
                  "(only for used BVLAN`s): ")
            print("-> {}".format(len(get_max_possible_path_a_to_b_filtered(
                input1, input2, set_bvls_in_use))))
        else:
            get_path_a_to_b(input1, input2)
            print("\nNumber of different pathes available " +
                  "(for all calculated BVLANs): ")
            print("-> {}".format(len(get_max_possible_path_a_to_b(
                input1, input2))))
    except Exception:
        key_error_node()


@decorator_cls_and_return_mainmenu
def menu4():
    input1 = input("Type number of pathes that should be available: ")
    try:
        print("List of connections with {} possibilities: ".format(input1))
        print("===========================================")
        if int(input1) not in d_ect_path.keys():
            print("There are no connections with {} possibilities.".format(
                  input1))
            print("Please check pathes summary.")
        else:
            for i in range(len(d_ect_path[int(input1)])):
                print(d_ect_path[int(input1)][i])
    except Exception:
        print("Please type a number!")
        input("Type anykey to continue")
        os.system("cls")


@decorator_cls_and_return_mainmenu
def menu10():
    print("Multicast states per node: ")
    print("==========================")
    for k in d_mcast_states:
        print("Node {:>{width}} has {:>5d} multicast states".format(
            k, len(d_mcast_states[k]), width=max_nodename))


@decorator_cls_and_return_mainmenu
def menu11():
    input1 = input("Select node: ")
    try:
        print("Multicast states of node {}".format(input1))
        x = PrettyTable()
        x.field_names = ["Source", "ISID", "BVLAN", "IIF", "OIL"]
        newlist = []  # To remove duplicates.
        for i in range(len(d_mcast_states[input1])):
            if d_mcast_states[input1][i] not in newlist:
                newlist.append(d_mcast_states[input1][i])
        for i in range(len(newlist)):
            x.add_row(newlist[i])
        x.sortby = "ISID"
        print(x)
    except Exception:
        key_error_node()


@decorator_cls_and_return_mainmenu
def menu20():
    print("List of all edges: ")
    print("=================  ")
    d_edges_weight_table = PrettyTable()
    d_edges_weight_table.field_names = ["Edge", "Metric"]
    for k, v in d_edges_weight.items():
        d_edges_weight_table.add_row([k, v])
    print(d_edges_weight_table)
    print("Total number of edges: {}".format(len(d_edges.keys())))


@decorator_cls_and_return_mainmenu
def menu21():
    d_edges_table = PrettyTable()
    d_edges_table.field_names = ["Edge", "Connections"]
    for key in d_edges.keys():
        d_edges_table.add_row([key, len(d_edges[key])])
    d_edges_table.sortby = "Connections"
    d_edges_table.reversesort = True
    print(d_edges_table)


@decorator_cls_and_return_mainmenu
def menu22():
    d_edges_connections = PrettyTable()
    d_edges_connections.field_names = ["Source", "Destination", "BVLAN"]
    input1 = input("Select source node: ")
    input2 = input("Select destination node: ")
    print("===========================================")
    edge = input1 + "___" + input2
    try:
        for i in d_edges[edge]:
            d_edges_connections.add_row([i[0], i[1], i[2]])
        print("Connections using edge {}".format(edge))
        print(d_edges_connections)

    except Exception:
        input3 = input("Edge does not exist! Select correct source and " +
                       "destination node. Display list of edges (y/n)?: ")
        if input3 == "y" or input3 == "Y":
            for key in d_edges.keys():
                print(key)
            input("\nType anykey to continue")
            os.system("cls")
        else:
            os.system("cls")


@decorator_cls_and_return_mainmenu
def menu30():
    apsp_ect_counter = 0
    print("Info: Every byte of SystemID that needs to be compared, \
will increment the value (i.e. one ECT decision with tiebreaker in 8 byte is \
counted as 8!")
    print("Tip : An earlier tiebreaker in SystemID will reduce the number of \
computations.\n")
    for k in d_node:
        print("Node {:>{width}} performed {:>4d} ECT computations \
(XOR comparisations) per BVLAN - hence in total: {}".format(
                 k, d_dijkstra[k][low_bvl_to_test].ectcnt,
                 d_dijkstra[k][low_bvl_to_test].ectcnt*len(set_bvls_in_use),
                 width=max_nodename))
        apsp_ect_counter += d_dijkstra[k][low_bvl_to_test].ectcnt*len(
            set_bvls_in_use)
    print("================================================")
    print("APSP would need {} ECT computations (XOR comparisations) ".format(
        apsp_ect_counter))
    # v2_0: Print d_tiebreaker_byte.
    print("\nTiebreaker statistics (position of SystemID it was found): ")
    print("Info: This is calculated for all BVLANs in range {} to {} \
(not only used BVLANs)!\n".format(low_bvl_to_test, max_bvl_to_test - 1))
    for k in d_tiebreaker_byte:
        print("Tiebreaker found in byte #{}: {}".format(
            k, d_tiebreaker_byte[k]))


@decorator_cls_and_return_mainmenu
def menu31():
    print("-> Time for all-pair-shortest-path (apsp) calculation: {}".format(
        time_apsp_calculation))
    print("-> Time for apsp calculation per BVLAN               : {}".format(
        time_apsp_calculation / (max_bvl_to_test - low_bvl_to_test)))
    print("-> Time for multicast-table population               : {}".format(
        time_mcast_calculation))
    print("\nInfo:")
    print("Multicast computation time might appear lower than experienced.")
    print("This is due to a remarkable time needed for output formatting.")


def menu_exit():
    os.system("exit")
    os.system("cls")


def main():

    # Assign dict with ECT to BVLAN.
    # (assuming BVLAN 4000 has ECT#1 and 4015 has ECT#15).
    global d_ect  # Global declaration of dict, as also accessed by other cls.
    d_ect = {"4000": b"\x00",
             "4001": b"\xff",
             "4002": b"\x88",
             "4003": b"\x77",
             "4004": b"\x44",
             "4005": b"\x33",
             "4006": b"\xcc",
             "4007": b"\xbb",
             "4008": b"\x22",
             "4009": b"\x11",
             "4010": b"\x66",
             "4011": b"\x55",
             "4012": b"\xaa",
             "4013": b"\x99",
             "4014": b"\xdd",
             "4015": b"\xee",
             }
    # Add experimental ECT ID´s, for calculating max possible path diversity.
    # Yeah right, VLAN´s > 4096 doesn´t exist - but that doesn´t matter here.
    for x in range(4016, 4256, 1):
        ectvalue = hex(x-4000)  # ECT value as string.
        ectvalue_hex = ectvalue.replace("0x", "")
        d_ect[str(x)] = bytes.fromhex(ectvalue_hex)

    global d_ect_path
    d_ect_path = {}
    global d_node
    global d_edges_weight
    d_edges_weight = {}
    global d_edges
    d_edges = {}
    global d_node_xor
    d_node_xor = {}
    # v2_0: Statistics, on which byte a tiebreaker was found for ECT.
    global d_tiebreaker_byte
    d_tiebreaker_byte = {0: 0, 1: 0, 2: 0, 3: 0, 4: 0, 5: 0, 6: 0, 7: 0}
    # To speed up computation in large networks, reduce here the number of
    # BVLANs to be tested (range is "low_bvl_to_test" - "num_bv_to_test" ).
    global max_bvl_to_test
    # max_bvl_to_test should be your highest BVLAN plus one.
    max_bvl_to_test = 4016
    global low_bvl_to_test
    low_bvl_to_test = 4000
    # Plausibility checks, for correct BVLAN range.
    assert max_bvl_to_test > low_bvl_to_test, "Please correct \
max_bvl_to_test and low_bvl_to_test. max value must be larger than low value!"
    assert 4256 >= max_bvl_to_test >= 4000, "Please correct \
max_bvl_to_test! Specify value between 4000 and 4256"
    assert 4256 >= low_bvl_to_test >= 4000, "Please correct \
low_bvl_to_test! Specify value between 4000 and 4256"

    # Plausibility check on d_node.
    # Check for unidirectional links and cost mismatch.
    for i in d_node:
        for j in d_node[i]["Links"]:
            unidirectional_link_check = True
            link_cost_mismatch = True
            for k in d_node[j[0]]["Links"]:
                if i == k[0]:
                    unidirectional_link_check = False
                if j[1] == k[1]:
                    link_cost_mismatch = False
            assert unidirectional_link_check is False,\
                "Unidirectional Link found on Node {}!".format(i)
            assert link_cost_mismatch is False,\
                "Link cost mismatch found on Node {}!".format(i)

    # Check for ISID/BVL mismatch and same BridgeID.
    for i in d_node:
        for isid in d_node[i]["ISID"]:
            check_isid = isid[0]
            check_bvl = isid[1]
            for j in d_node:
                if j != i:
                    assert d_node[i]["BridgeID"] != d_node[j]["BridgeID"],\
                     "Same BridgeID found on Node {} and {}!".format(i, j)
                for isid in d_node[j]["ISID"]:
                    if isid[0] == check_isid:
                        assert check_bvl == isid[1],\
                             "Wrong ISID/BVLAN relation on Node {}!".format(j)

    """ Modify your SPB network here to test behavior
    with more ISID`s / BVLAN`s etc. """

    # Temporary modification of d_node (i.e. to test more ISID`s):
    # for k in d_node.keys():
    #    for i in range(1000):
    #    d_node[k]["ISID"].append([100000+i, 4000+ i % 8, 0, 0])
    #    d_node[k]["ISID"].append([200000+i, 4000+i%256, 1, 1])10
    #    d_node[k]["ISID"].append([202020, 4015, 1, 1])
    #    d_node[k]["ISID"].append([202021, 4014, 1, 1])
    #    d_node[k]["ISID"].append([202022, 4013, 1, 1])
    #    d_node[k]["ISID"].append([202023, 4012, 1, 1])

    # Remove "off" nodes from dict.
    list_keys_to_remove = []
    for k in d_node:
        if d_node.get(k)["SystemState"] != "on":
            list_keys_to_remove.append(k)
    for k in range(len(list_keys_to_remove)):
        d_node.pop(list_keys_to_remove[k])

    # v2_0: Populate d_node_xor, results of xor operation of all BVLANs.
    for bvl in range(low_bvl_to_test, max_bvl_to_test):
        d_node_xor[bvl] = {}
        for k in d_node:
            d_node_xor[bvl][k] = []

    for bvl in range(low_bvl_to_test, max_bvl_to_test):
        for k in d_node:
            # Populate dict with xor result of SystemID.
            for i in range(8):
                d_node_xor[bvl][k].append([byte_xor(
                    d_node.get(k)["BridgeID"][i][0], d_ect.get(str(bvl)))])

    # Get maximum string lenght of node names.
    global max_nodename
    max_nodename = 0
    for k in d_node.keys():
        if len(k) > max_nodename:
            max_nodename = len(k)

    # Get a set of all BVLANs in loaded dict,
    # for having better display filter possibilites.
    global set_bvls_in_use
    set_bvls_in_use = set()
    for k in d_node:
        for i in range(len(d_node.get(k)["ISID"])):
            set_bvls_in_use.add(d_node.get(k)["ISID"][i][1])

    # Plausibility check, that all used BVLANs are calculated.
    for i in set_bvls_in_use:
        assert i in range(low_bvl_to_test, max_bvl_to_test), "Please correct \
max_bvl_to_test and low_bvl_to_test! BVLAN " + str(i) + " is used, bot not \
in that range!"

    # Create Graph and add all edges to d_edges and d_edges_weight.
    graph = Graph()
    for k in d_node:
        for i in range(len(d_node.get(k)["Links"])):
            if d_node.get(k)["SystemState"] == "on":
                # Only if node v is active (is a key in the d_node).
                a = d_node.get(k)["Links"][i][0]
                if a in d_node.keys():
                    graph.add_edge(k, d_node.get(k)["Links"][i][0],
                                   d_node.get(k)["Links"][i][1])
                    key_d_edges = str(k) + "___" + str(a)
                    d_edges[key_d_edges] = []
                    d_edges_weight[key_d_edges] = d_node.get(k)["Links"][i][1]

    global d_dijkstra  # Global declaration of dict, as accessed by other fct.
    d_dijkstra = {}  # Set every node as S.
    start_time_dijkstra = time.time()
    for n in d_node:
        d_dijkstra[n] = {}  # Initialize dict.
        for bvl in range(low_bvl_to_test, max_bvl_to_test):
            dijkstra = DijkstraSPF(graph, n, str(bvl))
            d_dijkstra[n][bvl] = dijkstra
    global time_apsp_calculation
    time_apsp_calculation = time.time() - start_time_dijkstra

    # Populate d_edges.
    for src, dst, bvl in product(d_node, d_node, set_bvls_in_use):
        for i in range(len(d_dijkstra[src][bvl].get_path(dst))):
            if len(d_dijkstra[src][bvl].get_path(dst)) > 0 and \
                  i < len(d_dijkstra[src][bvl].get_path(dst))-1:
                key_d_edges = d_dijkstra[src][bvl].get_path(dst)[i] + "___" +\
                              d_dijkstra[src][bvl].get_path(dst)[i+1]
                if key_d_edges in d_edges.keys():
                    d_edges[key_d_edges].append([src, dst, bvl])

    # Calculate tandem multicast mode.
    start_time_mcast = time.time()
    d_tandem_isid = {}
    set_bvl = set()
    d_tandem_mcast_src = {}
    d_tandem_mcast_rcv = {}
    global d_mcast_states
    d_mcast_states = {}  # End results to be displayed per node or in total.
    for node in d_node:
        d_mcast_states[node] = []

    # Populate d_tandem_mcast_src.
    for k in d_node:
        for i in range(len(d_node.get(k)["ISID"])):
            if d_node.get(k)["ISID"][i][2] == 1:  # Check T Bit of ISID.
                # Create the dict for MCSource and bvl.
                if k in d_tandem_mcast_src.keys():
                    d_tandem_mcast_src[k].\
                        append(d_node.get(k)["ISID"][i][0])
                else:
                    d_tandem_mcast_src[k] = \
                        [d_node.get(k)["ISID"][i][0]]
                if d_node.get(k)["ISID"][i][0] not in d_tandem_isid.keys():
                    d_tandem_isid[d_node.get(k)["ISID"][i][0]] = \
                                    [d_node.get(k)["ISID"][i][1]]
                    set_bvl.add(d_node.get(k)["ISID"][i][1])
                else:
                    pass  # Only add Tandem ISID to dict, if key not yet exist.

    # Populate d_tandem_mcast_rcv.
    for k in d_node:
        for i in range(len(d_node.get(k)["ISID"])):
            if d_node.get(k)["ISID"][i][3] == 1:  # Check R Bit of ISID.
                if k in d_tandem_mcast_rcv.keys():
                    d_tandem_mcast_rcv[k].append(d_node.get(k)["ISID"][i][0])
                else:
                    d_tandem_mcast_rcv[k] = \
                     [d_node.get(k)["ISID"][i][0]]
                if d_node.get(k)["ISID"][i][0] not in d_tandem_isid.keys():
                    d_tandem_isid[d_node.get(k)["ISID"][i][0]] = \
                     [d_node.get(k)["ISID"][i][1]]
                    set_bvl.add(d_node.get(k)["ISID"][i][1])
                else:
                    pass  # Only add Tandem ISID to dict, if key not yet exist.

    # Calculate MC States, for every node, in every BVLAN.
    for node, mc_src_a in product(d_node, d_tandem_mcast_src):
        for isid in d_tandem_mcast_src[mc_src_a]:  # For every isid.
            for bvl, mc_rcv_b in product(
             d_tandem_isid[isid], d_tandem_mcast_rcv):
                if isid in d_tandem_mcast_rcv[mc_rcv_b] and node in \
                 d_dijkstra[mc_src_a][bvl].get_path(mc_rcv_b)[1:-1]:
                    pos_node_path = d_dijkstra[mc_src_a][bvl].get_path(
                     mc_rcv_b).index(node)
                    iil = d_dijkstra[mc_src_a][bvl].\
                        get_path(mc_rcv_b)[pos_node_path-1]
                    oil = d_dijkstra[mc_src_a][bvl].\
                        get_path(mc_rcv_b)[pos_node_path+1]
                    d_mcast_states[node].append(
                        [mc_src_a, isid, bvl, iil, oil])
                elif isid in d_tandem_mcast_rcv[mc_rcv_b] and node in \
                    d_dijkstra[mc_src_a][bvl].get_path(mc_rcv_b)[0] and \
                        len(d_dijkstra[mc_src_a][bvl].get_path(mc_rcv_b)) > 1:
                    pos_node_path = d_dijkstra[mc_src_a][bvl].get_path(
                        mc_rcv_b).index(node)
                    iil = "-"  # Node is START of path.
                    oil = d_dijkstra[mc_src_a][bvl].\
                        get_path(mc_rcv_b)[pos_node_path+1]
                    d_mcast_states[node].append(
                     [mc_src_a, isid, bvl, iil, oil])
                elif isid in d_tandem_mcast_rcv[mc_rcv_b] and node in \
                    d_dijkstra[mc_src_a][bvl].get_path(
                     mc_rcv_b)[-1] and len(d_dijkstra[mc_src_a][bvl].get_path(
                         mc_rcv_b)) > 1:
                    pos_node_path = d_dijkstra[mc_src_a][bvl].get_path(
                        mc_rcv_b).index(node)
                    iil = d_dijkstra[mc_src_a][bvl].get_path(
                        mc_rcv_b)[pos_node_path-1]
                    oil = "-"  # Node is END of path.
                    d_mcast_states[node].append(
                        [mc_src_a, isid, bvl, iil, oil])
    global time_mcast_calculation
    time_mcast_calculation = time.time() - start_time_mcast

    # Consolidate OIL of d_mcast_states.
    # First step, remove entry`s with empty OIL, if:
    #    there are also entry´s with an OIL
    # That is the case, if node is receiver and also forkout point.

    for k in d_mcast_states.keys():
        mcast_states_to_delete = []  # Empty initialize for all nodes:
        for i in range(len(d_mcast_states[k])):
            if d_mcast_states[k][i][4] == "-":  # Check if node has emtpy OIL.
                # If so, check if node is also forkout point.
                for j in range(len(d_mcast_states[k])):
                    if d_mcast_states[k][i][0] == d_mcast_states[k][j][0] and \
                       d_mcast_states[k][i][1] == d_mcast_states[k][j][1] and \
                       d_mcast_states[k][i][2] == d_mcast_states[k][j][2] and \
                       d_mcast_states[k][i][3] == d_mcast_states[k][j][3] and \
                       d_mcast_states[k][i][4] != d_mcast_states[k][j][4]:
                        # Add entry with empty OIL to delete list.
                        mcast_states_to_delete.append(d_mcast_states[k][i])

        for x in range(len(mcast_states_to_delete)):
            if mcast_states_to_delete[x] in d_mcast_states[k]:
                d_mcast_states[k].remove(mcast_states_to_delete[x])

    # Second step, if OIL is not empty and has several entries - merge them.
    # That is the case, if node is forkout point with several OIL`s or sender.

    for k in d_mcast_states.keys():
        mcast_states_to_delete = []  # Empty initialize for all nodes:
        for i in range(len(d_mcast_states[k])):
            # Initialize the OIL that will be applied to:
            #     mcast state after iteration.
            new_oil = d_mcast_states[k][i][4]
            # Check that entry was not already checked.
            if d_mcast_states[k][i] not in mcast_states_to_delete:
                # Check if node has OIL (is a forkout point).
                if d_mcast_states[k][i][4] != "-":
                    for j in range(i+1, len(d_mcast_states[k])):
                        if d_mcast_states[k][i][0] == \
                            d_mcast_states[k][j][0] and \
                           d_mcast_states[k][i][1] == \
                            d_mcast_states[k][j][1] and \
                           d_mcast_states[k][i][2] == \
                            d_mcast_states[k][j][2] and \
                           d_mcast_states[k][i][3] == \
                           d_mcast_states[k][j][3]:
                            # Add duplicate entry to delete list.
                            mcast_states_to_delete.append(d_mcast_states[k][j])
                            # If OIL not already in OIL -> merge it.
                            if d_mcast_states[k][j][4] not in new_oil:
                                new_oil += ", " + d_mcast_states[k][j][4]
            # Finally, assign new OIL to mcast entry.
            d_mcast_states[k][i][4] = new_oil

        for x in range(len(mcast_states_to_delete)):
            if mcast_states_to_delete[x] in d_mcast_states[k]:
                d_mcast_states[k].remove(mcast_states_to_delete[x])

    #  Population of d_ect_path.
    for a in d_node:
        for b in d_node:
            k = len(get_max_possible_path_a_to_b(a, b))  # Get number of hops.
            if k in d_ect_path.keys():
                if (a + " -> " + b) not in d_ect_path[k] and a != b:
                    d_ect_path[k].append(a + " -> " + b)
                else:
                    pass  # Path already exist for this key.
            elif a != b:
                d_ect_path[k] = [a + " -> " + b]

    # Main User Menu for getting input from user.
    main_menu()


if __name__ == "__main__":
    main()
