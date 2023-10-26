import os
import sys

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

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"{sys.argv[0]} <chrom_id.lst> <chrAll.fna> <output_dir>")
        sys.exit(1)
    _, chrom_file, ref_file, out_dir = sys.argv


    if out_dir[-1] == "/":
        out_dir = out_dir[:-1]
    os.makedirs(out_dir, exist_ok=True)

    os.chdir(out_dir)

    chrom_ids = []
    with open(chrom_file, "r") as chrom_fd:
        for line in chrom_fd:
            chrom_id, left, right = line.strip().split("\t")
            chrom_ids.append((left, right, chrom_id))
        chrom_fd.close()
    
    ref_dict = get_fasta_todict(ref_file)

    afile=f"ref_a.fasta"
    bfile=f"ref_b.fasta"
    for (left, right, chrom_id) in chrom_ids:
        os.system("echo "" > " + afile)
        with open(afile, "w") as afd:
            afd.write(f">{left}\n")
            afd.write(ref_dict[left]+"\n")
            afd.close()
        
        os.system("echo "" > " + bfile)
        with open(bfile, "w") as bfd:
            bfd.write(f">{right}\n")
            bfd.write(ref_dict[right]+"\n")
            bfd.close()

        prefix=f"chr{chrom_id}"
        run_mummer1 = f"nucmer --mum --prefix={prefix} {afile} {bfile}"
        run_mummer2 = f"delta-filter -1 {prefix}.delta > {prefix}_1.delta"
        run_mummer3 = f"show-snps -x 10 -C -T -H -l -q {prefix}_1.delta > {prefix}.snp"

        os.system(run_mummer1)
        os.system(run_mummer2)
        os.system(run_mummer3)
    
    