import os
import sys

if __name__ == "__main__":
    if len(sys.argv) != 4:
        print(f"{sys.argv[0]} <gfa> <rename.map> <out_file>")
        sys.exit(1)
    
    _, gfa_file, map_file, out_file = sys.argv

    os.system("echo ""  > " + out_file)

    mapper = {}
    with open(map_file, "r") as mfd:
        for line in mfd:
            s, t = line.strip().split("\t")
            mapper[t] = s
        mfd.close()

    ofd = open(out_file, "w")
    with open(gfa_file, "r") as gfd:
        for line in gfd:
            splited = line.strip().split("\t")
            # only consider S or L to rename the ids
            if splited[0] == "S":
                splited[1] = mapper[splited[1]]
            elif splited[0] == "L":
                splited[1] = mapper[splited[1]]
                splited[3] = mapper[splited[3]]
            sstr = "\t".join(splited)
            ofd.write(f"{sstr}\n")
        gfd.close()
    ofd.close()

