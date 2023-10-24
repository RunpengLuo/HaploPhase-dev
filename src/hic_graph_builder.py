from lib.utils import *
from lib.file_io import *

def load_contracts(bed_file: str, matrix_file: str):
    num_entry = line_counter(bed_file)
    bed_arr = [None for _ in range(num_entry + 1)] # 1-based
    name_idx_dict = {}

    with open(bed_file, "r") as bed_fd:
        for line in bed_fd:
            name, start, end, idx = line.strip().split("\t")
            bed_arr[int(idx)] = (name, int(start), int(end))
            if name not in name_idx_dict:
                name_idx_dict[name] = []
            name_idx_dict[name].append(int(idx))
        bed_fd.close()

    name_st_dict = {}
    for name, arr in name_idx_dict.items():
        s, t = arr[0], arr[-1]
        m = (t - s + 1) // 2
        name_st_dict[name] = (s, s + m, t)
    inter_contacts = {}
    intra_contacts = {}
    for name in name_st_dict.keys():
        intra_contacts[name] = []
    with open(matrix_file, "r") as mat_fd:
        for line in mat_fd:
            uidx, vidx, val = line.strip().split("\t")
            link = (int(uidx), int(vidx), int(round(float(val))))
            uname = bed_arr[int(uidx)][0]
            vname = bed_arr[int(vidx)][0]
            if uname == vname:
                # intra contacts
                intra_contacts[uname].append(link)
            else:
                # inter contacts
                if (uname, vname) not in inter_contacts:
                    inter_contacts[(uname, vname)] = []
                inter_contacts[(uname, vname)].append(link)
        mat_fd.close()
    return bed_arr, name_st_dict, inter_contacts, intra_contacts

def build_hic_graph(clen_dict: dict, bed_file: str, matrix_file: str, out_dir: str):
    print("start building hic graph..")
    _, name_st_dict, inter_contacts, _ = load_contracts(bed_file, matrix_file)

    graph, vdict, edict = init_graph(False)
    for name in name_st_dict.keys():
        u = graph.add_vertex()
        graph.vp.id[u] = name
        vdict[name] = u
    
    for (uname, vname), links in inter_contacts.items():
        u = vdict[uname]
        v = vdict[vname]

        eval = sum([val for (_, _, val) in links])
        e = graph.add_edge(u, v)
        graph.ep.weight[e] = eval
        graph.ep.color[e] = "B"
        edict[(uname, vname)] = e

    out_file = out_dir + "/hic_graph.gfa"
    export_gfa_file(graph, vdict, edict, out_file, clen_dict)
    print("Done")
    return  out_file