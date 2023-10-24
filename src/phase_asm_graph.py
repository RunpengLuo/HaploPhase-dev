from graph_tool.all import Graph
from lib.file_io import *
from lib.utils import *
from lib.Steiner_Tree import *
from matching import get_subgraphs


"""
    build subgraph by subset vertices and edges.
    keep id and phase information.
"""
def build_subgraph(graph: Graph, vdict: dict, edict: dict, subsets: list):
    ograph, ovdict, oedict = init_graph(directed=False)
    ograph.vp.phase = ograph.new_vertex_property("int", val=0)

    for sid in subsets:
        u = ograph.add_vertex()
        ograph.vp.id[u] = sid
        ograph.vp.phase[u] = graph.vp.phase[vdict[sid]]
        ovdict[sid] = u
    for (uid, vid), old_e in edict.items():
        if uid in ovdict and vid in ovdict:
            new_e = ograph.add_edge(ovdict[uid], ovdict[vid])
            ograph.ep.color[new_e] = "B" # default color
            oedict[(uid, vid)] = new_e

    return ograph, ovdict, oedict

def print_phase(graph: Graph, vdict: dict, p: int, prefix=""):
    s = []
    for uid, u in vdict.items():
        if graph.vp.phase[u] == p:
            s.append(uid)
    
    print(prefix+f"phase {p}: ", ",".join(s))

# propogate phasing label along simple path
def label_propogation(graph: Graph, vdict: dict, edict: dict):
    is_simple = {}
    for uid, u in vdict.items():
        degree = u.out_degree()
        if (uid, uid) in edict: # ignore self-loop
            degree -= 1
        is_simple[uid] = degree <= 2

    num_propogation = -1
    while num_propogation != 0:
        num_propogation = 0
        for (uid, vid), e in edict.items():
            u, v = e.source(), e.target()
            if graph.vp.phase[u] == -1:
                if is_simple[uid] and graph.vp.phase[v] != -1:
                    graph.vp.phase[u] = graph.vp.phase[v]
                    num_propogation += 1
            
            if graph.vp.phase[v] == -1:
                if is_simple[vid] and graph.vp.phase[u] != -1:
                    graph.vp.phase[v] = graph.vp.phase[u]
                    num_propogation += 1

        print("number of propogation: ", num_propogation)
    return

def phase_asm(contig_dict: dict, gfa_file: str, hic_graph_file: str, bin0file: str, bin1file: str, out_dir: str):
    
    bin0dict = get_fasta_todict(bin0file)
    bin1dict = get_fasta_todict(bin1file)

    hic_graph, hic_vdict, hic_edict = import_gfa_file(hic_graph_file)
    graph, vdict, edict = import_gfa_file(gfa_file, ignore_weight=True) # color = B
    graph.vp.phase = graph.new_vertex_property("int", val=-1)
    # assign phase label to vertices
    for vid, v in vdict.items():
        if vid in bin0dict:
            graph.vp.phase[v] = 0
        elif vid in bin1dict:
            graph.vp.phase[v] = 1
        else:
            graph.vp.phase[v] = -1

    node_weight0 = dict.fromkeys(vdict.keys(), 0) # hic connection to phase 0
    node_weight1 = dict.fromkeys(vdict.keys(), 0)
    for uid in vdict.keys():
        for p0id in bin0dict.keys():
            if (uid, p0id) in hic_edict:
                node_weight0[uid] += hic_graph.ep.weight[hic_edict[(uid, p0id)]]
            elif (p0id, uid) in hic_edict:
                node_weight0[uid] += hic_graph.ep.weight[hic_edict[(p0id, uid)]]
            else:
                pass
        for p1id in bin1dict.keys():
            if (uid, p1id) in hic_edict:
                node_weight1[uid] += hic_graph.ep.weight[hic_edict[(uid, p1id)]]
            elif (p1id, uid) in hic_edict:
                node_weight1[uid] += hic_graph.ep.weight[hic_edict[(p1id, uid)]]
            else:
                pass
        print(f"Node {uid} - {node_weight0[uid]} - {node_weight1[uid]}")
    
    edge_weight0 = dict.fromkeys(edict.keys(), 0)
    edge_weight1 = dict.fromkeys(edict.keys(), 0)
    for (uid, vid) in edict.keys():
        # be careful about the reverse sign
        edge_weight0[(uid, vid)] = node_weight1[uid] + node_weight1[vid]
        edge_weight1[(uid, vid)] = node_weight0[uid] + node_weight0[vid]

    print_phase(graph, vdict, 0, "Before label propagation")
    print_phase(graph, vdict, 1, "Before label propagation")

    label_propogation(graph, vdict, edict)

    print_phase(graph, vdict, 0, "After label propagation")
    print_phase(graph, vdict, 1, "After label propagation")
    
    graph.vp.group = graph.new_vertex_property("int", val=0)
    gid_max = get_subgraphs(graph, vdict, edict)

    phase0 = []
    phase1 = []

    for gid in range(1, gid_max):
        gnames = [graph.vp.id[ver] for ver in graph.vertices() if graph.vp.group[ver] == gid]
        print(f"{gid}th subgraph: ", ",".join(gnames))

        phase_mask0 = [uid for uid in gnames if graph.vp.phase[vdict[uid]] != 0]
        phase_mask1 = [uid for uid in gnames if graph.vp.phase[vdict[uid]] != 1]

        # build a subgraph with no phase 0
        sgraph0, svdict0, sedict0 = build_subgraph(graph, vdict, edict, phase_mask1)
        sgraph0.vp.group = sgraph0.new_vertex_property("int", val=0)
        s0gid_max = get_subgraphs(sgraph0, svdict0, sedict0)
        for s0gid in range(1, s0gid_max):
            s0gnames = [sgraph0.vp.id[ver] for ver in sgraph0.vertices() if sgraph0.vp.group[ver] == s0gid]
            ssgraph0, ssvdict0, ssedict0 = build_subgraph(sgraph0, svdict0, sedict0, s0gnames)

            steiner_tree0 = steiner_tree(ssgraph0, ssvdict0, ssedict0, contig_dict, 0, edge_weight0)
            if steiner_tree0 != None:
                for uid in steiner_tree0:
                    if graph.vp.phase[vdict[uid]] == -1:
                        # undefined node assigned to at least one phase via steiner tree
                        graph.vp.phase[vdict[uid]] = 2 # mark involved
                    phase0.append(uid)
        
        # build a subgraph with no phase 0
        sgraph1, svdict1, sedict1 = build_subgraph(graph, vdict, edict, phase_mask0)
        sgraph1.vp.group = sgraph1.new_vertex_property("int", val=0)
        s1gid_max = get_subgraphs(sgraph1, svdict1, sedict1)
        for s1gid in range(1, s1gid_max):
            s1gnames = [sgraph1.vp.id[ver] for ver in sgraph1.vertices() if sgraph1.vp.group[ver] == s1gid]
            ssgraph1, ssvdict1, ssedict1 = build_subgraph(sgraph1, svdict1, sedict1, s1gnames)

            steiner_tree1 = steiner_tree(ssgraph1, ssvdict1, ssedict1, contig_dict, 1, edge_weight1)
            if steiner_tree1 != None:
                for uid in steiner_tree1:
                    if graph.vp.phase[vdict[uid]] == -1:
                        # undefined node assigned to at least one phase via steiner tree
                        graph.vp.phase[vdict[uid]] = 2 # mark involved
                    phase1.append(uid)

    print("[After steiner tree] Phase 0: ", ",".join(phase0))
    print("[After steiner tree] Phase 1: ", ",".join(phase1))

    # TODO those are coming from the subgraphs with no terminals or not traversed by st path
    # manually add isolated contig to the closest bin if available
    # cov_cutoff = 99
    # depth_cutoff = 15 # simply homozygous coverage // 4, move to argv later
    i_value_delta = 0

    for cid, ctg in vdict.items():
        if graph.vp.phase[ctg] != -1:
            continue
        i_value_a = node_weight0[cid]
        i_value_b = node_weight1[cid]
        if max(i_value_a, i_value_b) <= i_value_delta:
                continue

        if abs(i_value_a - i_value_b) <= i_value_delta:
            continue
        if i_value_a > i_value_b:
            # goto a
            phase0.append(cid)
            graph.vp.phase[vdict[cid]] = 0
            print(f"add isolated contig {cid} to phase 0 with weight {i_value_a}")
        else:
            # goto b
            phase1.append(cid)
            graph.vp.phase[vdict[cid]] = 1
            print(f"add isolated contig {cid} to phase 1 with weight {i_value_b}")

    print("[Final] Phase 0: ", ",".join(phase0))
    print("[Final] Phase 1: ", ",".join(phase1))

    out_bin0 = out_dir + "/tree_bin0.fasta"
    out_bin1 = out_dir + "/tree_bin1.fasta"
    out_binN = out_dir + "/tree_binN.fasta"
    Create(out_bin0)
    Create(out_bin1)
    Create(out_binN)

    fd0 = open(out_bin0, "w")
    fd1 = open(out_bin1, "w")
    fdN = open(out_binN, "w")
    for id0 in phase0:
        fd0.write(f">{id0}\n{contig_dict[id0]}\n")
    for id1 in phase1:
        fd1.write(f">{id1}\n{contig_dict[id1]}\n")
    for uid, u in vdict.items():
        if graph.vp.phase[u] == -1:
            fdN.write(f">{uid}\n{contig_dict[uid]}\n")
    fd0.close()
    fd1.close()
    fdN.close()

    return out_bin0, out_bin1, out_binN
