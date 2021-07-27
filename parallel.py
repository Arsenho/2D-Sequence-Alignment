import utils
import mpi4py.MPI as MPI
import time
import argparse
import sys
import numpy as np

comm = MPI.COMM_WORLD
nprocs = comm.Get_size()
rank = comm.Get_rank()

parser = argparse.ArgumentParser(sys.argv[0])
parser.add_argument('--rows', type=int, nargs=1, help="Number of rows")

parser.add_argument('--columns', type=int, nargs=1, help="Number of columns")

args = parser.parse_args()

if args.rows is None:
    sys.stderr.write("Missing minimum k\n")
    sys.exit(1)

if args.columns is None:
    sys.stderr.write("Missing Number of equivalent classes per processors\n")
    sys.exit(1)

pre_treatment_time = 0.0


def t_matrice_func(x, y):
    start_pre_treatment = time.time()

    rows_motif_1 = len(x)
    rows_motif_2 = len(y)

    columns_motif_1 = len(x[0])
    columns_motif_2 = len(y[0])

    my_motif_1 = utils.split_data(x, nprocs, rank)
    # print(rank, " my motif length = ", my_motif_1)

    comm.barrier()
    dr_mat = utils.dr_matrice_func(x)
    dc_mat = utils.dc_matrice_func(x)
    ir_mat = utils.ir_matrice_func(y)
    ic_mat = utils.ic_matrice_func(y)

    r_mat = utils.r_matrice_func(my_motif_1, y)
    c_mat = utils.c_matrice_func(my_motif_1, y)
    root = 0
    row = 0

    comm.barrier()

    # if rank != 0:
    #     row = len(my_motif_1)
    #     comm.send(row, dest=root, tag=rank)
    #     comm.send(r_mat, dest=root, tag=rank)
    #     comm.send(c_mat, dest=root, tag=rank)
    # else:
    #     all_r_mat = []
    #     all_c_mat = []
    #     # 4 dimension matrix initializations
    #     for i in range(rows_motif_1):
    #         inter_1 = []
    #         for j in range(columns_motif_1):
    #             inter_2 = []
    #             for k in range(rows_motif_2):
    #                 inter_3 = []
    #                 for l in range(columns_motif_2):
    #                     inter_3.append(0)
    #                 inter_2.append(inter_3)
    #             inter_1.append(inter_2)
    #         all_r_mat.append(inter_1)
    #         all_c_mat.append(inter_1)
    #
    #     for rk in range(1, nprocs):
    #         row = comm.recv(source=rk, tag=rk)
    #
    #         r_mat_prime = comm.recv(source=rk, tag=rk)
    #         # print(rank, " ", r_mat_prime[:1])
    #         for i in range(row):
    #             for j in range(columns_motif_1):
    #                 for k in range(rows_motif_2):
    #                     for l in range(columns_motif_2):
    #                         all_r_mat[i][j][k][l] = r_mat_prime[i][j][k][l]
    #
    #         c_mat_prime = comm.recv(source=rk, tag=rk)
    #         for i in range(row):
    #             for j in range(columns_motif_1):
    #                 for k in range(rows_motif_2):
    #                     for l in range(columns_motif_2):
    #                         all_c_mat[i][j][k][l] = c_mat_prime[i][j][k][l]

    # del c_mat
    # del r_mat

    comm.barrier()
    pre_treatment_time = time.time() - start_pre_treatment

    # if rank == 0:
    start_treatment = time.time()
    rows_my_motif_1 = len(my_motif_1)

    t_matrix = np.zeros((rows_motif_1, columns_motif_1, rows_motif_2, columns_motif_2), dtype=int)

    # Initalization maginales of T table
    for i in range(rows_my_motif_1):
        for j in range(columns_motif_1):
            for k in range(rows_motif_2):
                for l in range(columns_motif_2):
                    t_matrix[i][j][k][l] = 0
                    t_matrix[0][j][k][l] = (k + 1) * (l + 1)
                    t_matrix[i][0][k][l] = (k + 1) * (l + 1)
                    t_matrix[i][j][0][l] = (i + 1) * (j + 1)
                    t_matrix[i][j][k][0] = (i + 1) * (j + 1)


    # Filling of T table
    for i in range(1, rows_my_motif_1):
        for j in range(1, columns_motif_1):
            for k in range(1, rows_motif_2):
                for l in range(1, columns_motif_2):
                    t_matrix[i][j][k - 1][l] = max(
                        t_matrix[i - 1][j][k][l] + dr_mat[i][rows_motif_1 - 1],
                        t_matrix[i][j - 1][k][l] + dc_mat[i][rows_motif_1 - 1],
                        t_matrix[i][j][k - 1][l] + ir_mat[i][rows_motif_2 - 1],
                        t_matrix[i][j][k - 1][l - 1] + ic_mat[i][rows_motif_2 - 1],
                        t_matrix[i - 1][j][k - 1][l] + r_mat[i][j][k][l],
                        t_matrix[i][j - 1][k][l - 1] + c_mat[i][j][k][l],
                        t_matrix[i - 1][j - 1][k - 1][l - 1] + c_mat[i - 1][j][k - 1][l] + r_mat[i][j][k][l],
                        t_matrix[i - 1][j - 1][k - 1][l - 1] + c_mat[i][j][k][l] + r_mat[i][j - 1][k][l - 1]
                    )

    # for i in range(rows_motif_1):
    #     for j in range(columns_motif_1):
    #         for k in range(rows_motif_2):
    #             for l in range(columns_motif_2):
    #                 print("[{}][{}][{}][{}] = ".format(i, j, k, l), t_matrice[i][j][k][l])

    treatment_time = time.time() - start_treatment

    comm.barrier()
    global_pre_treatment_time = comm.reduce(pre_treatment_time, op=MPI.SUM, root=0)

    if rank == 0:
        return global_pre_treatment_time / nprocs, treatment_time
    else:
        return pre_treatment_time


rows = args.rows[0]
columns = args.columns[0]

motif_1 = utils.get_motif(rows, columns)
motif_2 = utils.get_motif(rows, columns)

start = time.time()
result = t_matrice_func(motif_1, motif_2)
end = time.time() - start

if rank == 0:
    dtf = [result[0], result[1], end, rows, columns, nprocs]
    stats_file_name = "stats/parallel"

    stat_file_name = 'stats/parallel/stats_for_parallel_rows_{}_columns_{}.csv'.format(rows, columns)
    utils.list_to_same_file(stat_file_name, [dtf])

    print(rank, " Computations = ", result)
    print(rank, " Global computation time = ", end)
