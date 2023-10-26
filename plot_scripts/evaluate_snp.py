import os
import sys
from transformer import *

import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

def parse_mfile(mfile: str, k=21):
    """mummer snp detection results - DEBUG only

    Args:
        mfile (str): _description_
    """

    indels = 0
    non_isolated_snp = 0
    higher_snp = 0
    snp = 0

    ignoreN = 0

    A_snps = set()
    B_snps = set()
    snp_index = k//2
    with open(mfile, "r") as fd:
        for line in fd:
            splits = line[:-1].split("\t")
            [sub1, sub2] = splits[1:3]
            [kmer1, kmer2] = splits[8:10]
            if not is_alpha(sub1) or not is_alpha(sub2):
                indels += 1
                continue
            if 'N' in kmer1 or 'N' in kmer2:
                ignoreN += 1
                continue

            A_snps.add(kmer1)
            A_snps.add(reverse_alpha(kmer1))
            B_snps.add(kmer2)
            B_snps.add(reverse_alpha(kmer2))


            if alpha2int(kmer1[:snp_index] + kmer1[snp_index+1:]) != alpha2int(kmer2[:snp_index] + kmer2[snp_index+1:]):
                if hamming_distance(kmer1, kmer2) <= 2:
                    non_isolated_snp += 1
                else:
                    higher_snp += 1
            else:
                snp += 1
    
        fd.close()

    print("total: ", indels + non_isolated_snp + higher_snp + snp)
    print("SNP: ", snp)
    print("non_isolated SNP: ", non_isolated_snp)
    print("higher order SNP: ", higher_snp)
    print("indels: ", indels)
    print("ignore N's kmer: ", ignoreN)

    return A_snps, B_snps

def get_fasta_todict(fasta_file: str):
    res = {}
    with open(fasta_file, "r") as fa_fd:
        sid = ""
        seq = ""
        for line in fa_fd:
            if line.startswith(">"):
                # process previous entry
                if sid != "":
                    res[sid] = seq
                sid = line.strip()[1:].split()[0]
                seq = ""
            else:
                seq += line.strip()
        fa_fd.close()
        if sid != "":
            res[sid] = seq
    return res

def get_snp_count(hap_file: str, Asnp_global: set, Bsnp_global: set):

    hapAdict = get_fasta_todict(hap_file)
    hapA_xB = []
    hapA_yA = []
    counter = 0
    for cid, cseq in hapAdict.items():
        counter += 1
        if counter % 50 == 0:
            print("finished: ", counter)
        yAsnp = 0
        xBsnp = 0
        len_ctg_seq = len(cseq)
        for sub_i in range(len_ctg_seq - k + 1):
            sub = cseq[sub_i : sub_i + k]
            # check global snp 
            if sub in Asnp_global:
                yAsnp += 1
            elif sub in Bsnp_global:
                xBsnp += 1
            else:
                pass
        if yAsnp != 0 or xBsnp != 0:
            hapA_xB.append(max(xBsnp//1000, 0))
            hapA_yA.append(max(yAsnp//1000, 0))
    return hapA_xB, hapA_yA

if __name__ == "__main__":
    if len(sys.argv) != 8:
        print(f"{sys.argv[0]} <chrom_id.lst> <chrAll.fna> <mummer_dir> <k> <HapA.fasta> <HapB.fasta> <prefix>")
        sys.exit(1)
    _, chrom_file, chrAll_file, mum_dir, k, hapA_file, hapB_file, prefix = sys.argv

    k = int(k)

    if mum_dir[-1] == "/":
        out_dir = mum_dir[:-1]

    if prefix[-1] == "_":
        prefix = prefix[:-1]

    # get minimap2 result
    if not os.path.exists(prefix+"a2ref_A.paf"):
        exec_mmp = f"minimap2 -x asm5 -t 16 --secondary=no {chrAll_file} {hapA_file} > {prefix}+a2ref_A.paf"
        os.system(exec_mmp)
    
    if not os.path.exists(prefix+"a2ref_B.paf"):
        exec_mmp = f"minimap2 -x asm5 -t 16 --secondary=no {chrAll_file} {hapB_file} > {prefix}+a2ref_B.paf"
        os.system(exec_mmp)

    chr_id2name = {}
    name2chr_id = dict()
    with open(chrom_file, "r") as chrom_fd:
        for line in chrom_fd:
            chrom_id, left, right = line.strip().split("\t")
            chr_id2name[chrom_id] = (left, right)
            name2chr_id[left] = chrom_fd
            name2chr_id[right] = chrom_fd
        chrom_fd.close()
    
    # get SNPs per chromosome
    Asnp_global = set()
    Asnps = []
    Bsnp_global = set()
    Bsnps = []
    for chr_id, (left, right) in chr_id2name.items():
        print(f"Mummer hash table for chr{chr_id}")
        Asnp, Bsnp = parse_mfile(f"{mum_dir}/chr{chr_id}.snp")
        print(len(Asnp), len(Bsnp))
        Asnp_global = Asnp_global.union(Asnp)
        Bsnp_global  = Bsnp_global.union(Bsnp)

        itsect = Asnp.intersection(Bsnp)
        Asnp = Asnp.difference(itsect)
        Bsnp = Bsnp.difference(itsect)

        Asnps.append(Asnp)
        Bsnps.append(Bsnp)

    # mutual exclusive
    itsect = Asnp_global.intersection(Bsnp_global)
    Asnp_global = Asnp_global.difference(itsect)
    Bsnp_global = Bsnp_global.difference(itsect)

    print(len(itsect))
    print(len(Asnp_global))
    print(len(Bsnp_global))

    # process Haplotype A's contigs
    # hapA_xB, hapA_yA = get_snp_count(hapA_file, Asnp_global, Bsnp_global)
    hapA_xB, hapA_yA = get_snp_count(hapA_file, Asnps[0], Bsnps[0])
    print(hapA_xB)
    print(hapA_yA)

    # hapB_xB, hapB_yA = get_snp_count(hapB_file, Asnp_global, Bsnp_global)
    hapB_xB, hapB_yA = get_snp_count(hapB_file, Asnps[0], Bsnps[0])
    print(hapB_xB)
    print(hapB_yA)

    # df = pd.DataFrame({"# Haplotype A k-mers": hapA_yA, "# Haplotype B k-mers": hapA_xB})
    fig, ax =plt.subplots(1,1)
    plt.scatter(hapA_xB, hapA_yA, label="hap1", c="red")
    plt.scatter(hapB_xB, hapB_yA, label="hap2", c="blue")
    
    _max1 = max(max(hapA_xB), max(hapA_yA))
    _max2 = max(max(hapB_xB), max(hapB_yA))
    _max = max(_max1, _max2)

    plt.yticks(np.arange(0, _max, 20))
    plt.xticks(np.arange(0, _max, 20))
    ax.set_xlabel("# Haplotype B k-mers")
    ax.set_ylabel("# Haplotype A k-mers")
    ax.set_title("Hap-mer Plot")
    x0,x1 = ax.get_xlim()
    y0,y1 = ax.get_ylim()
    ax.set_aspect(abs(x1-x0)/abs(y1-y0))
    ax.legend()
    # plt.show()
    plt.savefig(prefix + "_hapmer_plot.png")

    # hapB_xB, hapB_yA = get_snp_count(hapB_file, Asnp_global, Bsnp_global)

    
