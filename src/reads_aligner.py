import os
import sys
from lib.utils import *
from lib.file_io import *

import pysam

def classify_reads_samtools(contig_dict: dict, bam_file: str, out_dir: str):
    print("classify reads to contig..")
    rid_dir = out_dir+"/reads"
    os.makedirs(rid_dir)
    for rname in contig_dict.keys():
        rname_sep = rname.split()[0]
        lst_file = rid_dir + f"/sep_{rname_sep}.lst"
        # UNMAP,SECONDARY,QCFAIL,DUP
        System(f"samtools view -F 0x704 {bam_file} {rname_sep} | cut -f1 | sort | uniq > {lst_file}")
    print("done")
    return

def classify_reads(contig_dict: dict, bam_file: str, out_dir: str):
    print("classify reads to contig..")
    samfile = pysam.AlignmentFile(bam_file, "rb")

    rid_sdict = {}
    print("#rname\tnum_entries\treads_fail\treads_pass\n")
    for rname in contig_dict.keys():
        pass_qnames = set()
        fail_qnames = set()

        rname_sep = rname.split()[0]
        iter = list(samfile.fetch(rname_sep))
        for aln in iter:
            # UNMAP,SECONDARY,QCFAIL,DUP
            if aln.is_unmapped or aln.is_secondary or aln.is_qcfail or aln.is_duplicate:
                fail_qnames.add(aln.query_name)
                continue
            pass_qnames.add(aln.query_name)
        rid_sdict[rname] = list(pass_qnames)
        print(f"{rname_sep}\t{len(iter)}\t{len(fail_qnames)}\t{len(pass_qnames)}")

    rid_dir = out_dir+"/reads"
    os.makedirs(rid_dir)

    for cid in contig_dict.keys():
        cid_sep = cid.split()[0]
        # create name.lst first
        name_lst_file = rid_dir + f"/sep_{cid_sep}.lst"
        Create(name_lst_file)
        with open(name_lst_file, "w") as lst_fd:
            for rid in rid_sdict[cid]:
                lst_fd.write(f"{rid}\n")
            lst_fd.close()
        # rid_file = rid_dir + f"/sep_{cid_sep}.fastq"
        # System(f"seqtk subseq {read_file} {name_lst_file} > {rid_file}")

        # print(f"{cid_sep} processed, stored {len(rid_sdict[cid])} reads")
    print("Done")
    return rid_dir

def concat_lsts(bin_fa: str, bin_lst: str, rid_dir: str):
    bin_set = set()
    bin_dict = get_edges(bin_fa)
    for cid in bin_dict.keys():
        cid_sep = cid.split()[0]
        lst_file = rid_dir + f"/sep_{cid_sep}.lst"
        assert (os.path.exists(lst_file))
        s = set()
        with open(lst_file, "r") as lst_fd:
            for line in lst_fd:
                s.add(line.strip())
        lst_fd.close()
        bin_set = bin_set.union(s)
    
    with open(bin_lst, "w") as lst_fd:
        for rid in bin_set:
            lst_fd.write(f"{rid}\n")
        lst_fd.close()
    return bin_lst, bin_set

def bin_reads(read_file: str, bin0_fa: str, bin1_fa: str, binN_fa: str, out_dir: str):
    print("binning reads...")

    rid_dir = out_dir + "/reads"

    bin0_lst, bin0_set = concat_lsts(bin0_fa, out_dir + "/bin0_rnames.lst", rid_dir)
    bin0_fq = out_dir + "/bin0_reads.fastq"
    System(f"seqtk subseq {read_file} {bin0_lst} > {bin0_fq}")

    bin1_lst, bin1_set = concat_lsts(bin1_fa, out_dir + "/bin1_rnames.lst", rid_dir)
    bin1_fq = out_dir + "/bin1_reads.fastq"
    System(f"seqtk subseq {read_file} {bin1_lst} > {bin1_fq}")

    binN_lst, binN_set = concat_lsts(binN_fa, out_dir + "/binN_rnames.lst", rid_dir)
    binN_fq = out_dir + "/binN_reads.fastq"
    System(f"seqtk subseq {read_file} {binN_lst} > {binN_fq}")

    binAll = bin0_set.union(bin1_set).union(binN_set)
    print("bin 0: ", len(bin0_set))
    print("bin 1: ", len(bin1_set))
    print("bin N: ", len(binN_set))
    print("bin All: ", len(binAll))

    print("Done")
    return bin0_fq, bin1_fq


def rm_dup_reads(bin_fq: str, out_dir: str):

    return
