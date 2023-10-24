from lib.utils import *
from lib.file_io import *

def build_double_graph(clen_dict: dict, snp_matching_file: str, hic_graph_file, out_dir: str):
    # vset: conflict contig pairs from snp_matching_file
    conf_pairs = list(from_matching_file(snp_matching_file).keys())
    # eset: edges from hic graph exclude edges within the pairs

    hic_graph, hic_vdict, hic_edict = import_gfa_file(hic_graph_file)
    bi_graph, bi_vdict, bi_edict = init_graph(False)
    for (uid, vid) in conf_pairs:
        # only consider the conflict edge connected vertices
        u = bi_graph.add_vertex()
        bi_graph.vp.id[u] = uid
        bi_vdict[uid] = u
        v = bi_graph.add_vertex()
        bi_graph.vp.id[v] = vid
        bi_vdict[vid] = v
    for (uid, vid), e in hic_edict.items():
        if uid in bi_vdict and vid in bi_vdict:
            if min_pair(uid, vid) in conf_pairs:
                continue # SNP overlay edge, ignore
            bi_e = bi_graph.add_edge(bi_vdict[uid], bi_vdict[vid])
            bi_graph.ep.weight[bi_e] = hic_graph.ep.weight[e]
            bi_graph.ep.color[bi_e] = "B"
            bi_edict[(uid, vid)] = bi_e
    
    out_file = out_dir + "/bi_graph.gfa"
    export_gfa_file(bi_graph, bi_vdict, bi_edict, out_file, clen_dict)
    print("Done")

    return out_file