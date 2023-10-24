from lib.utils import *
from lib.file_io import *
from lib.transformer import *

def build_snp_graph(contig_dict: dict, clen_dict: dict, snp_file: str, snp_ext_file: str, out_dir: str):
    print("start building snp graph..")
    snp0, snp1, iso_cov, k = get_snp_file(snp_file) # TODO
    num_snps = len(snp0)
    print("number of snps: ", num_snps)
    snp_groups0, snp_groups1, _ = get_snp_extend_file(snp_ext_file, k)


    kmer_dict = dict()
    # kmer_dict[k_i] = [G_0, G_1, ...] indicates k_i is co-exists in groups G_0, G_1,...
    # 0/1 for isolated SNP L/R, 2/3 for non-isolated SNP L/R
    for idx, groups in enumerate([snp_groups0, snp_groups1]):
        for i, group in enumerate(groups):
            for j, kmer in enumerate(group):
                if kmer not in kmer_dict:
                    kmer_dict[kmer] = []
                kmer_dict[kmer].append((idx, i, j))
    print("kmer dict constructed")

    uniq_kmer_dict = dict()
    for km, idds in kmer_dict.items():
        if len(idds) == 1:
            alpha=int2alpha(km,k)
            uniq_kmer_dict[alpha] = idds[0]
            uniq_kmer_dict[reverse_alpha] = idds[0]



    snp_ctg_arr_x = [[] for _ in range(num_snps)]
    snp_ctg_arr_y = [[] for _ in range(num_snps)]
    snp_ctg_arr = [snp_ctg_arr_x, snp_ctg_arr_y]
    # check all contigs
    print(f"performing snp-contig perfect matching..")
    for ctg_id, ctg_seq in contig_dict.items():
        len_ctg_seq = len(ctg_seq)
        for sub_i in range(len_ctg_seq - k + 1):
            sub = ctg_seq[sub_i : sub_i + k]
            if sub in uniq_kmer_dict:
                (xy_idx, i_idx, _) = uniq_kmer_dict[sub]
                snp_ctg_arr[xy_idx][i_idx].append(ctg_id)
    
    
    conflict_dict = {}
    intermediate_file = out_dir + "/snp2ctg.csv"
    Create(intermediate_file)
    with open(intermediate_file, "w") as out_fd:
        for i in range(num_snps):
            fst_terms, snd_terms = set(snp_ctg_arr_x[i]), set(snp_ctg_arr_y[i])
            l_str = ",".join(fst_terms)
            r_str = ",".join(snd_terms)
            out_fd.write(f"{l_str},-,{r_str}\n")

            if len(fst_terms) == 1 and len(snd_terms) == 1:
                # current SNP uniquely align to two 
                # (potentially same) contigs on each part
                lctg_id, rctg_id = list(fst_terms)[0], list(snd_terms)[0]
                if lctg_id != rctg_id:
                    # two distinct contigs
                    mp = min_pair(lctg_id, rctg_id)
                    conflict_dict[mp] = conflict_dict.get(mp, 0) + 1
        out_fd.close()
    
    intermediate_file2 = out_dir + "/snp_conflict_edges.csv"
    Create(intermediate_file2)
    with open(intermediate_file2, "w") as conf_fd:
        for (uid, vid), val in conflict_dict.items():
            conf_fd.write(f"{uid},{vid},{val}\n")
        conf_fd.close()

    # construct the snp graph
    graph, vdict, edict = init_graph(False)
    for uid in clen_dict.keys():
        unode = graph.add_vertex()
        graph.vp.id[unode] = uid
        vdict[uid] = unode

    for (uid, vid), w in conflict_dict.items():
        assert (uid, vid) not in edict
        e = graph.add_edge(vdict[uid], vdict[vid])
        graph.ep.weight[e] = w
        graph.ep.color[e] = "B"
        edict[(uid, vid)] = e

    out_file = out_dir + "/snp_graph.gfa"
    export_gfa_file(graph, vdict, edict, out_file, clen_dict)
    print("Done")
    return out_file