import os
import sys

def System(command: str):
    return os.system(command)

def Create(file: str):
    return System("echo "" > " + file)

def line_counter(file: str):
    c = 0
    with open(file, "r") as fd:
        for _ in fd:
            c += 1
        fd.close()
    return c

def get_edges(fasta_file):
    """
    assume sequence is in one-line
    """
    unitigs = {}
    prev_id = ""
    with open(fasta_file, "r") as fd:
        for l in fd:
            if l.startswith(">"):
                prev_id = l.strip().split()[0][1:]
            else:
                unitigs[prev_id] = l.strip()
        fd.close()

    return unitigs

def get_file(file: str, is_fasta: bool):
    if is_fasta:
        return get_fasta(file)
    else:
        return get_fastq(file)

def get_fasta(fasta_file: str):
    res = []
    with open(fasta_file, "r") as fa_fd:
        sid = ""
        seq = ""
        for line in fa_fd:
            if line.startswith(">"):
                # process previous entry
                if sid != "":
                    res.append((sid, seq))
                sid = line.strip()[1:]
                seq = ""
            else:
                seq += line.strip()
        fa_fd.close()
        if sid != "":
            res.append((sid, seq))
    return res

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
                sid = line.strip()[1:]
                seq = ""
            else:
                seq += line.strip()
        fa_fd.close()
        if sid != "":
            res[sid] = seq
    return res

def get_fastq(fastq_file: str):
    reads = []
    counter = 0
    with open(fastq_file, "r") as fd:
        prev_id = ""
        prev_seq = ""
        i = 0
        for l in fd:
            if i == 0:
                assert l.startswith("@")
                counter += 1

                if prev_id != "":
                    reads.append((prev_id, prev_seq))
                
                # read id, get first non-space id
                prev_id = l.strip()
                prev_seq = ""
            elif i == 1:
                prev_seq += l.strip()
            elif i == 2:
                if not l == "+\n":
                    # multiline sed
                    prev_seq += l.strip()
                    continue
            i = (i+1) % 4
        fd.close()
        if prev_id != "":
            reads.append((prev_id, prev_seq))

    print("total reads: ", counter)
    return reads

def reformat_fasta(in_file: str, out_file: str):
    os.system("echo "" > " + out_file)
    with open(in_file, "r") as in_fd:
        with open(out_file, "w") as out_fd:
            sid = ""
            seq = ""
            for line in in_fd:
                if line.startswith(">"):
                    # process previous entry
                    if sid != "":
                        for i in range(len(seq)):
                            if seq[i] not in {"A","C","T","G","N"}:
                                print(sid)
                                print("invalid: ", seq[i:i+5])
                        out_fd.write(sid + "\n")
                        out_fd.write(seq + "\n")
                    sid = line.strip()
                    
                    seq = ""
                else:
                    seq += line.strip().upper()
            if sid != "":
                out_fd.write(sid + "\n")
                out_fd.write(seq + "\n")
            out_fd.close()
        in_fd.close()

def store_reads(reads: dict, only_id=True):
    # store basic read informations
    System("echo "" > reads.txt")
    with open("reads.txt", "w") as fd:
        for rid, rseq in reads.items():
            if only_id:
                fd.write(f"{rid}\n")
            else:
                fd.write(f"{rid}:{rseq}\n")
        fd.close()

def get_reads(only_id=True):
    ids = []
    with open("reads.txt", "r") as fd:
        for line in fd:
            if only_id:
                ids.append(line.strip())
            else:
                raise Exception
        fd.close()
    return ids

def get_ctuple(cid: str):
    t, idx = cid.split("_")
    return (int(t[:-1]), t[-1], int(idx))
