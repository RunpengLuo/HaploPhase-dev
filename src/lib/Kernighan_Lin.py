from graph_tool.all import Graph

def compute_d_values(graph: Graph, vdict: dict, edict: dict, s_bisection: dict):
    # compute d values
    d_values = {}
    for vid, vnode in vdict.items():
        # compute d value
        my_phase = s_bisection[vid]
        i_value = 0
        e_value = 0
        for e in vnode.all_edges():
            neighbor = e.source() if e.target() == vnode else e.target()
            if s_bisection[graph.vp.id[neighbor]] != my_phase:
                # external connection
                e_value += graph.ep.weight[e]
            else:
                # internal connection
                i_value += graph.ep.weight[e]
        d_values[vid] = e_value - i_value
    return d_values

def exchange(set_a: list, set_b: list, s_bisection: dict, max_k: int, visited_k: list, num_flip: int, g_values: list):
    # perform exchange on first max_k operations
    for i in range(max_k):
        curr_a, curr_b = set_a[visited_k[i]], set_b[visited_k[i]]
        a_phase, b_phase = s_bisection[curr_a], s_bisection[curr_b]
        s_bisection[curr_b] = a_phase
        s_bisection[curr_a] = b_phase
        print(f"{curr_a}\t{a_phase}->{b_phase}\t{curr_b}\t{b_phase}->{a_phase}\t{g_values[i]}")
    for i in range(num_flip):
        a, b = set_a[i], set_b[i]
        if s_bisection[a] != 0:
            # perform a flip, useless though, just for visualisation
            set_a[i] = b
            set_b[i] = a
    return

def compute_cost(graph: Graph, edict: dict, set_a: list, set_b: list):
    left_cost = 0
    right_cost = 0
    external_cost = 0
    for (uid, vid), e in edict.items():
        if (uid in set_a and vid in set_b) or (vid in set_a and uid in set_b):
            external_cost += graph.ep.weight[e] 
            # print(f"{uid}-{vid}\t{graph.ep.weight[e]}")
        if (uid in set_a and vid in set_a):
            left_cost += graph.ep.weight[e]
        if (uid in set_b and vid in set_b):
            right_cost += graph.ep.weight[e]
    return left_cost, external_cost, right_cost  


def get_weight(graph: Graph, edict: dict, uid, vid):
    if (uid, vid) in edict:
        return graph.ep.weight[edict[(uid, vid)]]
    elif (vid, uid) in edict:
        return graph.ep.weight[edict[(vid, uid)]]
    else:
        return 0

"""
edge marked as color "C" is conflicting edge
"""
def kernighan_lin(graph: Graph, vdict: dict, edict: dict, s_bisection: dict, set_a: list, set_b: list):
    print("Before Greedy: ", compute_cost(graph, edict, set_a, set_b))
    num_flip = len(set_a)
    num_iteration = 100
    while num_iteration != 0:
        num_iteration -= 1
        # compute d values
        d_values = compute_d_values(graph, vdict, edict, s_bisection)
        
        g_values = []
        visited_k = []

        for _ in range(num_flip):
            # pick a flip that maximize g = D[a] + D[b] - 2w(a,b) for a in A and b in B
            max_g_value = 0
            max_k = -1
            for k in range(num_flip):
                if k not in visited_k:
                    a = set_a[k]
                    b = set_b[k]
                    curr_g = d_values[a] + d_values[b] - 2 * get_weight(graph, edict, a, b)
                    # w(a,b) should be 0 in here
                    if curr_g > max_g_value:
                        max_k = k
                        max_g_value = curr_g
            if max_k == -1: # TODO add a minimum swap threshold
                print("no more g-value reduction at this point, skip")
                break
            max_a = set_a[max_k]
            max_b = set_b[max_k]

            # remove a and b from further consideration in this pass
            visited_k.append(max_k)
            g_values.append(max_g_value)

            a_phase = s_bisection[max_a]
            # update d_values, simulate flip but not truly flip
            for sid, s_phase in s_bisection.items():
                if sid == max_a or sid == max_b:
                    continue
                if s_phase == a_phase:
                    d_values[sid] = d_values[sid] + 2 * get_weight(graph, edict, sid, max_a) - 2 * get_weight(graph, edict, sid, max_b)
                else:
                    # b_phase
                    d_values[sid] = d_values[sid] + 2 * get_weight(graph, edict, sid, max_b) - 2 * get_weight(graph, edict, sid, max_a)

        g_max = 0
        max_k = 0
        for i in range(1, num_flip):
            acc_val = sum(g_values[:i])
            if acc_val > g_max:
                g_max = acc_val
                max_k = i
        
        print(f"Maximum g_value: {g_max} by {max_k}th iteration")
        if g_max > 0:
            exchange(set_a, set_b, s_bisection, max_k, visited_k, num_flip, g_values)
            print("Greedy: ", compute_cost(graph, edict, set_a, set_b))
        else:
            print("Final Greedy: ", compute_cost(graph, edict, set_a, set_b))
            # d-value should both be negative,
            return d_values

def optimal_bisection(graph: Graph, vdict: dict, edict: dict, set_a: list, set_b: list):
    gt_a = []
    gt_b = []

    for i in range(len(set_a)):
        a, b = set_a[i], set_b[i]
        a_phase = a.split("_")[0][-1]
        b_phase = b.split("_")[0][-1]
        if a_phase == "A":
            gt_a.append(a)
        else:
            gt_b.append(a)
        if b_phase == "A":
            gt_a.append(b)
        else:
            gt_b.append(b)

    # compute external cost
    print("multi Optimal: ", compute_cost(graph, edict, gt_a, gt_b))