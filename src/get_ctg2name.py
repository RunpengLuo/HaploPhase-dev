import os
import sys

from src.lib.utils import get_edges

if __name__ == "__main__":
    if len(sys.argv) != 7:
        print(f"{sys.argv[0]} <contig_file> <chrom_id.lst> <gtruth> <out_ctg2name>")
        sys.exit(1)
    _, contig_file, chrom_file, gtruth_file, out_file = sys.argv

    chrom_ids = {}
    with open(chrom_file, "r") as chrom_fd:
        for line in chrom_fd:
            chrom_id, left, right = line.strip().split("\t")
            chrom_ids[left] = [chrom_id + "A", 0]
            chrom_ids[right] = [chrom_id + "B", 0]
        chrom_fd.close()

    ctg2ref = {}
    with open(gtruth_file, "r") as gt_fd:
        for row in gt_fd:
            # each row is haplotype-related contigs
            rname, cnames = row.strip().split(": ")
            cname_arr = cnames.split(",")
            for cname in cname_arr:
                ctg2ref[cname] = rname
        gt_fd.close()
    
    contig_dict = get_edges(contig_file)
    ctg2name = {}
    # use ref_dict and chrom_ids to rename contig ids
    for name in contig_dict.keys():
        rname = ctg2ref[name]
        new_name, new_idx = chrom_ids[rname]
        chrom_ids[rname] = [new_name, new_idx + 1]
        ctg2name[name] = f"{new_name}_{new_idx}"
    
    with open(out_file, "w") as cfd:
        for name, new_name in ctg2name.items():
            cfd.write(f"{name}\t{new_name}\n")
        cfd.close()