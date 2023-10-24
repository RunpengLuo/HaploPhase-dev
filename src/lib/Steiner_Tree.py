import sys
from graph_tool.all import Graph

def path_weight(path: list, edge_weight: dict):
    pw = 0
    for i in range(len(path) - 1):
        uid = path[i]
        vid = path[i+1]
        if (uid, vid) in edge_weight:
            pw += edge_weight[(uid, vid)]
        else:
            pw += edge_weight[(vid, uid)]

    return len(path)

def steiner_tree(graph: Graph, vdict: dict, edict: dict, contig_dict: dict, phase: int, edge_weight: dict):
    print("Build Steiner tree for phase: ", phase)
    phases = [uid for uid, u in vdict.items() if graph.vp.phase[u] == phase]
    if len(phases) == 0:
        print("skip, no terminals be found")
        return None

    terminals = dict.fromkeys(phases, False)
    st_path_dict = {}
    # SV^2
    for uid in terminals.keys():
        # compute shortest path between u to all, Dijkstra
        # FIXME edge weight maybe assigned as the i-value.
        dist = dict.fromkeys(vdict.keys(), sys.maxsize)
        prev = dict.fromkeys(vdict.keys(), None)
        dist[uid] = 0

        Q = dict.fromkeys(vdict.keys(), True)
        while len(Q) != 0:
            vid,_ = min(Q.items(), key=lambda t: dist[t[0]])
            Q.pop(vid)

            for w in set(vdict[vid].out_neighbors()):
                wid = graph.vp.id[w]
                if wid in Q:
                    # alt = dist[vid] + 1 # w(v,w) = 1 for now
                    edge_w = 0
                    if (vid, wid) in edge_weight:
                        edge_w = edge_weight[(vid, wid)]
                    else:
                        edge_w = edge_weight[(wid, vid)]
                    alt = dist[vid] + edge_w
                    if alt < dist[wid]:
                        dist[wid] = alt
                        prev[wid] = vid
        
        # extract st path for each pair
        st_paths = {}
        for target in vdict.keys():
            st_paths[target] = []
            vid = target
            if prev[vid] != None:
                while vid != None:
                    st_paths[target].insert(0, vid)
                    vid = prev[vid]
            else:
                # self
                st_paths[target].insert(0, vid)
        st_path_dict[uid] = st_paths

    num_terminals = len(terminals)
    vset = set() # store id

    # init steiner tree with one terminal vertex, select the longest one.
    root = max(terminals, key=lambda k: len(contig_dict[k]))
    vset.add(root)
    terminals[root] = True
    num_terminals -= 1

    while num_terminals != 0:
        # select a closest terminal to the current tree
        min_st_path = []
        min_terminal = None
        min_weight = sys.maxsize
        for terminal, added in terminals.items():
            if added:
                continue
            # check shortest path between terminal and nodes in vset
            for vid in vset:
                path = st_path_dict[terminal][vid]
                pw = path_weight(path, edge_weight)
                if pw < min_weight:
                    min_st_path = path
                    min_terminal = terminal
                    min_weight = pw
        
        num_terminals -= 1
        terminals[min_terminal] = True
        # add to tree
        for vid in min_st_path:
            vset.add(vid)

        print(f"adding {min_terminal}-path with weight {min_weight}: ", ",".join(min_st_path))
        print()
    
    print("Steiner  tree: ", ",".join(vset))
    return list(vset)