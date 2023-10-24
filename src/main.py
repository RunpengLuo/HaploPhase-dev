import os
import sys
import time

from lib.utils import *

import hic_graph_builder as hic_graph_builder
import snp_graph_builder as snp_graph_builder
import double_graph_builder as double_graph_builder
import matching as matching
import greedy_bisection  as gbisection
import phase_asm_graph


if __name__ == "__main__":
    if len(sys.argv) != 8:
        print(f"{sys.argv[0]} <gfa_file> <asm_file> <bed> <matrix> <snp_file> <snp_ext_file> <output_dir>")
        sys.exit(1)
    _, gfa_file, asm_file, bed_file, matrix_file, snp_file, snp_ext_file, out_dir = sys.argv

    ts = time.time_ns()

    if out_dir[-1] == "/":
        out_dir = out_dir[:-1]
    os.makedirs(out_dir, exist_ok=True)

    # get all the contig informations, seq, id, length
    contig_dict = get_edges(asm_file)

    c_size_file = out_dir+"/contig_size.csv"
    clen_dict = {}
    Create(c_size_file)
    with open(c_size_file, "w") as cfd:
        for cid, cseq in contig_dict.items():
            cfd.write(f"{cid},{len(cseq)}\n")
            clen_dict[cid] = len(cseq)
        cfd.close()

    # snp file
    # snp_ext_file
    # -> snp graph
    ts_snp_graph=time.time_ns()
    if not os.path.exists(out_dir + "/snp_graph.gfa"):
        snp_graph_file = snp_graph_builder.build_snp_graph(contig_dict, clen_dict, snp_file, snp_ext_file, out_dir)
    else:
        print("snp graph file exists")
        snp_graph_file = out_dir + "/snp_graph.gfa"
    te_snp_graph=time.time_ns()
    print(f"[SNP graph construction] wall time elapsed: {round((te_snp_graph - ts_snp_graph) / (10 ** 9), 3)}s")

    # -> conflict contig pair (maximum weighted matching)
    ts_snp_matching=time.time_ns()
    if not os.path.exists(out_dir + "/matching.csv"):
        snp_matching_file = matching.max_weighted_matching(snp_graph_file, out_dir)
    else:
        print("snp graph file exists")
        snp_matching_file = out_dir + "/matching.csv"
    te_snp_matching=time.time_ns()
    print(f"[SNP matching] wall time elapsed: {round((te_snp_matching - ts_snp_matching) / (10 ** 9), 3)}s")

    # bed
    # matrix
    # -> hic graph
    ts_hic_graph=time.time_ns()
    if not os.path.exists(out_dir + "/hic_graph.gfa"):
        hic_graph_file = hic_graph_builder.build_hic_graph(clen_dict, bed_file, matrix_file, out_dir)
    else:
        print("hic graph file exists")
        hic_graph_file = out_dir + "/hic_graph.gfa"

    te_hic_graph=time.time_ns()
    print(f"[HiC graph construction] wall time elapsed: {round((te_hic_graph - ts_hic_graph) / (10 ** 9), 3)}s")

    # hic graph + conflict contig pair overlay
    # -> bisect graph
    ts_db_graph=time.time_ns()
    if not os.path.exists(out_dir + "/bi_graph.gfa"):
        double_graph_file = double_graph_builder.build_double_graph(clen_dict, snp_matching_file, hic_graph_file, out_dir)
    else:
        print("bi graph file exists")
        double_graph_file = out_dir + "/bi_graph.gfa"
    te_db_graph=time.time_ns()
    print(f"[DB graph construction] wall time elapsed: {round((te_db_graph - ts_db_graph) / (10 ** 9), 3)}s")
    
    # -> bisection
    ts_bisection=time.time_ns()
    if not os.path.exists(out_dir + "/bin0.fasta") or not os.path.exists(out_dir + "/bin1.fasta") or not os.path.exists(out_dir + "/binN.fasta"):
        bin0_file, bin1_file, binN_file = gbisection.bisection(contig_dict, clen_dict, double_graph_file, snp_matching_file, hic_graph_file, out_dir)
    else:
        print("bisection result exists")
        bin0_file = out_dir + "/bin0.fasta"
        bin1_file = out_dir + "/bin1.fasta"
        binN_file = out_dir + "/binN.fasta"
    te_bisection=time.time_ns()
    print(f"[Bisection] wall time elapsed: {round((te_bisection - ts_bisection) / (10 ** 9), 3)}s")
    
    # asm graph + hic graph + bin 0 + bin 1
    # -> final bisection via Steiner tree
    ts_steiner=time.time_ns()
    if not os.path.exists(out_dir + "/tree_bin0.fasta") or not os.path.exists(out_dir + "/tree_bin1.fasta") or not os.path.exists(out_dir + "/tree_binN.fasta"):
        tbin0_file, tbin1_file, tbinN_file = phase_asm_graph.phase_asm(contig_dict, gfa_file, hic_graph_file, bin0_file, bin1_file, out_dir)
    else:
        print("tree phasing exists")
        tbin0_file = out_dir + "/tree_bin0.fasta"
        tbin1_file = out_dir + "/tree_bin1.fasta"
        tbinN_file = out_dir + "/tree_binN.fasta"
    te_steiner=time.time_ns()
    print(f"[Steiner tree] wall time elapsed: {round((te_steiner - ts_steiner) / (10 ** 9), 3)}s")

    te=time.time_ns()
    print(f"total wall time elapsed: {round((te - ts) / (10 ** 9), 3)}s")
    sys.exit(0)

