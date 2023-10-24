from graph_tool.all import Graph
from lib.utils import *
from lib.transformer import *


def init_graph(directed=False):
    graph = Graph(directed)
    vdict = {}
    edict = {}
    graph.vp.id = graph.new_vertex_property("string", val="")
    graph.ep.weight = graph.new_edge_property("int", val=0)
    graph.ep.color = graph.new_edge_property("string", val="B")  
    return graph, vdict, edict

def export_gfa_file(graph: Graph, vdict: dict, edict: dict, filename: str, clen_dict: dict):
    # convert to GFA file
    Create(filename)
    with open(filename, "w") as gfa_fd:
        for vname in vdict.keys():
            gfa_fd.write(
                "S\t{0}\t{1}\tLN:i:{2}\n".format(
                    vname, "*", clen_dict[vname]
                )
            )
        for (uid, vid), e in edict.items():
            if graph.ep.color[e] != "B":
                continue
            gfa_fd.write(
            "L\t{0}\t{1}\t{2}\t{3}\t{4}M\n".format(
                uid, "+", vid, "+", graph.ep.weight[e])
            )
        gfa_fd.close()
    return

def import_gfa_file(gfa_file: str, ignore_weight=False):
    graph, vdict, edict = init_graph(False)

    with open(gfa_file, "r") as gfa_fd:
        for line in gfa_fd.readlines():
            if line.startswith("S"):
                # segment
                _, name = line.strip().split()[:2]
                vtx = graph.add_vertex()
                graph.vp.id[vtx] = name
                vdict[name] = vtx 
            elif line.startswith("L"):
                # link
                splits = line.strip().split()
                if len(splits) == 6:
                    [_, uname, _, vname, _, val] = splits
                else:
                    [_, uname, _, vname] = splits[:4]
                    assert ignore_weight # doesn't support advanced graph feature yet.

                e = graph.add_edge(vdict[uname], vdict[vname])
                if not ignore_weight:
                    graph.ep.weight[e] = int(val[:-1])
                else:
                    graph.ep.weight[e] = 0
                graph.ep.color[e] = "B"
                edict[(uname,vname)] = e
            else:
                # print("unsupported graph type", line[0])
                pass
        gfa_fd.close()
    return graph, vdict, edict

def get_snp_file(snp_file: str):
    k = 0
    iso_l_snps = []
    iso_r_snps = []
    iso_snp_covs = []
    with open(snp_file, "r") as snp_fd:
        for line in snp_fd.readlines():
            if line == "\n":
                break
            snp1, snp2, cov = line[:-1].split(" ")
            iso_l_snps.append(canonical_alpha2int(snp1))
            iso_r_snps.append(canonical_alpha2int(snp2))
            iso_snp_covs.append(int(cov))
            if k == 0:
                k = len(snp1)
        snp_fd.close()
    # print("Finish snp file parsing, Total SNP pair: ", len(iso_l_snps))
    return iso_l_snps, iso_r_snps, iso_snp_covs, k

def get_snp_extend_file(snp_ext_file: str, k):
    """Split extend SNP pair to fixed kmers
    """
    kmer_l_groups = []
    kmer_r_groups = []
    niso_snp_covs = []
    with open(snp_ext_file, "r") as ext_fd:
        i = 0
        for line in ext_fd.readlines():
            if line == "\n":
                break
            ext1, ext2, cov = line.strip().split(" ")
            kmer_l_groups.append([])
            kmer_r_groups.append([])
            niso_snp_covs.append(int(cov))
            for sub_i in range(len(ext1) - k + 1):
                sub1 = ext1[sub_i : sub_i + k]
                sub2 = ext2[sub_i : sub_i + k]
                if sub1.count("N") != 0 or sub2.count("N") != 0:
                    continue

                kmer_l_groups[i].append(canonical_alpha2int(sub1))
                kmer_r_groups[i].append(canonical_alpha2int(sub2))
            i += 1
    print("Finish snp extend file parsing")
    return kmer_l_groups, kmer_r_groups, niso_snp_covs

def from_matching_file(matching_file: str):
    mdict = {}
    with open(matching_file, "r") as fd:
        for line in list(fd.readlines())[1:]:
            uid, vid, w = line.strip().split(",")
            mdict[(uid, vid)] = int(w)
        fd.close()
    return mdict