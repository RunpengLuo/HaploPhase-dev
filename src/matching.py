from graph_tool.all import Graph
from graph_tool.topology import max_cardinality_matching

from lib.file_io import *
from lib.utils import *

def get_subgraphs(graph: Graph, vdict: dict, edict: dict):
    """
    assign group id for vertices, only using black edge
    """
    def bfs_coloring(start_name: str, group_idx: int):
        queue = [vdict[start_name]]
        while queue != []:
            curr = queue.pop()
            graph.vp.group[curr] = group_idx
            for adj_e in curr.all_edges():
                if graph.ep.color[adj_e] != "B":
                    continue
                next = adj_e.source() if adj_e.target() == curr else adj_e.target()
                if graph.vp.group[next] == 0 and next not in queue: # unassigned
                    queue.append(next)
        return

    gidx = 1
    for uname, u in vdict.items():
        if graph.vp.group[u] == 0:
            bfs_coloring(uname, gidx)
            gidx += 1
    print("groups: ", gidx)
    return gidx

def max_weighted_matching(graph_file: str, out_dir: str):
    print("maximum weighted matching..")
    graph, vdict, edict = import_gfa_file(graph_file)
    
    edge_matching = max_cardinality_matching(graph, weight=graph.ep.weight, edges=True)
    # remove non matching edges by marking the color as Grey
    for (uname, vname) in edict.keys():
        if not edge_matching[edict[(uname, vname)]]:
            graph.ep.color[edict[(uname, vname)]] = "G"
    graph.vp.group = graph.new_vertex_property("int", val=0)
    gid_max = get_subgraphs(graph, vdict, edict)

    matching_dict = {}
    for gid in range(1, gid_max):
        gnames = [graph.vp.id[ver] for ver in graph.vertices() if graph.vp.group[ver] == gid]
        if len(gnames) != 2:
            assert len(gnames) == 1
        else:
            [uid, vid] = gnames
            mp = min_pair(uid, vid)
            w = graph.ep.weight[edict[mp]]
            matching_dict[mp] = w
    
    out_file = out_dir + "/matching.csv"
    Create(out_file)
    with open(out_file, "w") as fd:
        fd.write(f"UID,VID,WEIGHT\n")
        for (uid, vid), w in matching_dict.items():
            fd.write(f"{uid},{vid},{w}\n")
        fd.close()

    print("Done")
    return out_file
