from lib.utils import *
from lib.file_io import *
from lib.Kernighan_Lin import *

def bisection(contig_dict: dict, clen_dict: dict, double_graph_file: str, snp_matching_file: str, hic_graph_file: str, out_dir: str):
    conf_dict = from_matching_file(snp_matching_file)
    conf_pairs = list(conf_dict.keys())
    graph, vdict, edict = import_gfa_file(double_graph_file)
    hic_graph, _, hic_edict = import_gfa_file(hic_graph_file)

    set_a = [a for (a,_) in conf_pairs]
    set_b = [b for (_,b) in conf_pairs]

    # build a random starting bisection
    # TODO optimize initial solution
    s_bisection = dict.fromkeys(vdict.keys(), -1)
    for (a,b) in conf_pairs:
        s_bisection[a] = 1
        s_bisection[b] = 0
        # s_bisection[a] = random.randint(0,1)
        # s_bisection[b] = (s_bisection[a] + 1) % 2

    # get bisection
    d_values = kernighan_lin(graph, vdict, edict, s_bisection, set_a, set_b)

    # manual fix FP conflicting pairs by checking d-values
    # improve the bisection by re-assign positive d-value contig if-any
    pruning(s_bisection, d_values)

    # scale d_values by contig length
    s_values = {}
    for cname, dval in d_values.items():
        s_values[cname] = round(dval / clen_dict[cname], 5)

    # display the bisection results
    out_stat_bi = out_dir + "/bi_stat.csv"
    Create(out_stat_bi)
    with open(out_stat_bi, "w") as osfd:
        osfd.write("uname,u-excluded,vname,v-excluded,ulen,vlen,g-value,d-value(u),s-value(u),d-value(v),s-value(v),HiC,SNP\n")
        for (a,b) in zip(set_a, set_b):
            g_value = d_values[a] + d_values[b] # w(a.b) should be 0
            hic_w = get_weight(hic_graph, hic_edict, a, b)
            snp_w = conf_dict[min_pair(a,b)]
            ulen = clen_dict[a]
            vlen = clen_dict[b]
            osfd.write(f"{a},{s_bisection[a] == 2},{b},{s_bisection[b] == 2},{ulen},{vlen},{g_value},{d_values[a]},{s_values[a]},{d_values[b]},{s_values[b]},{hic_w},{snp_w}\n")
        osfd.close() 

    for cid in contig_dict.keys():
        if cid not in s_bisection:
            s_bisection[cid] = 2 # TBD

    out_bin0 = out_dir + "/bin0.fasta"
    out_bin1 = out_dir + "/bin1.fasta"
    out_binN = out_dir + "/binN.fasta"
    Create(out_bin0)
    Create(out_bin1)
    Create(out_binN)

    fd0 = open(out_bin0, "w")
    fd1 = open(out_bin1, "w")
    fdN = open(out_binN, "w")
    for vid, dec in s_bisection.items():
        {0: fd0, 1: fd1, 2: fdN}[dec].write(f">{vid}\n{contig_dict[vid]}\n")
    fd0.close()
    fd1.close()
    fdN.close()

    return out_bin0, out_bin1, out_binN

def pruning(s_bisection: dict, d_values: dict):
    d_value_delta = 10 # minimal d_value cutoff, notice that d_value-> -inf    
    # flip positive d-value nodes to opponent
    # remove nodes lower than cutoff

    print(f"D-value delta: [{d_value_delta}, {d_value_delta}]")
    # re-assign according to d-values
    # -inf ~ -d_value_delta ~ d_value_delta ~ inf
    for uid in s_bisection.keys():
        d_value = d_values[uid]
        if d_value >= -d_value_delta and d_value <= d_value_delta:
            # not strongly favor to either nucleus bin, mark as contigs that deal later
            s_bisection[uid] = 2
            print(f"Manual Ignore: {uid},{d_value}")
        elif d_value > d_value_delta:
            # either no effect or should be placed to alternative bin but was restricted by conflict edge
            prev_assign = s_bisection[uid]
            new_assign = (prev_assign + 1) % 2
            s_bisection[uid] = new_assign
            print(f"Manual Fix: {uid},{d_value}: {prev_assign}->{new_assign}")
    return