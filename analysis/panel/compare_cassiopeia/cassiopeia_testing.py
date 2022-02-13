import networkx as nx
from networkx.algorithms.shortest_paths.dense import reconstruct_path
import cassiopeia as cas
import numpy as np
import pandas as pd
import time
import argparse
import faulthandler
import os

## Authored by Wexiang Fang
## revised for updated cassiopeia version by Abel Sapirstein

def reonstructTree(data_name, output_name, hybrid, replicate):
    df = pd.read_csv(data_name, sep = "\t", header = None)
    df = df.set_index(df.shape[1]-1)
    df = df.applymap(lambda x: int(x.strip('R-')))
    cas_tree = cas.data.CassiopeiaTree(character_matrix=df)
    hybrid.solve(cas_tree, logfile=f'logs/{replicate}_hybrid.log')
    vanilla_greedy.solve(cas_tree, collapse_mutationless_edges=True)
    rep=cas.data.to_newick(cas_tree.get_tree_topology())
    with open(output_name, "w") as f:
        f.write(rep)

def recon_tree(data_name, output_name):
    df = pd.read_csv(data_name, sep = "\t", header = None)
    df = df.set_index(df.shape[1]-1)
    df = df.applymap(lambda x: int(x.strip('R-')))
    l = [Node.Node(index, row.tolist()) for index, row in df.iterrows()]
    tree = ls.solve_lineage_instance(l, method="hybrid", threads = 1)
    pp.assign_samples_to_charstrings(tree[0].network, df)
    nw = convert_network_to_newick_format(tree[0].network, use_intermediate_names = False)
    with open(output_name, "w") as f:
        f.write(nw)

if __name__ == '__main__':
    faulthandler.enable()
    parser = argparse.ArgumentParser(description='Test for Expected SPR Behavior, currently configured to look at SPR and Node Sliding. \
        It does this by iterspersing bouts of node sliding with long MCMC chains of SPR moves')
    parser.add_argument('--treeNumber',"-t", type=int, required= True,help='which tree number to run - 1-12')
    parser.add_argument('--hgrnaNum',"-hg", type=int, required= True,help='25,50,100')
    args = parser.parse_args()
    data_dir = "../../6state_exp_test_data_v2/data/"
    max_index = args.treeNumber
    hgNum= args.hgrnaNum
    out_dir = "../outputs/cassiopeia_longer"
    run_time_file = "../outputs/cassiopeia_longer/6state_exp_test_cassi_run_time.csv"
    if not os.path.exists(os.path.dirname(run_time_file)):
        try:
            os.makedirs(os.path.dirname(run_time_file))
        except OSError as exc: # Guard against rare condition
            raise
    print("Creating Solvers")
    ilp_solver = cas.solver.ILPSolver()
    vanilla_greedy = cas.solver.VanillaGreedySolver()
    hybrid_solver = cas.solver.HybridSolver(top_solver=vanilla_greedy, bottom_solver=ilp_solver, cell_cutoff=40)
    print("Created Solvers .... opening runtime files")
    run_time_file = open("../outputs/cassiopeia_longer/6state_exp_test_cassi_run_time_2.txt", "a")
    run_time_array = np.array([])
    i = max_index
    i = str(i+1)
    string_num = i.zfill(4)
    t_time = time.time()
    reonstructTree(data_name = data_dir + i.zfill(4)+"_"+str(hgNum) + "_barcodes.txt",
        output_name = out_dir + i.zfill(4) +"_"+str(hgNum) + "_recontructed.txt",hybrid=hybrid_solver, replicate=string_num)
    run_time = time.time() - t_time
    run_time_file.write(f"{i},{hgNum},{str(run_time)}\n")
    run_time_array = np.append(run_time_array, run_time)

    run_time_file.close()
    np.savetxt(f"../outputs/cassiopeia_longer/6state_exp_test_cassi_run_time_bak_50_{i}.txt", run_time_array)
