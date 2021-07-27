import utils
import time
import argparse
import sys
import numpy as np

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


# Filling of T table
def t_matrice_func(x, y):
    rows_motif_1 = len(x)
    rows_motif_2 = len(y)

    columns_motif_1 = len(x[0])
    columns_motif_2 = len(y[0])

    start_pre_treatment = time.time()
    dr_mat = utils.dr_matrice_func(x)
    dc_mat = utils.dc_matrice_func(x)
    ir_mat = utils.ir_matrice_func(y)
    ic_mat = utils.ic_matrice_func(y)

    r_mat = utils.r_matrice_func(x, y)
    c_mat = utils.c_matrice_func(x, y)
    pre_treatment_time = time.time() - start_pre_treatment

    start_treatment = time.time()
    t_matrix = np.zeros((rows_motif_1, columns_motif_1, rows_motif_2, columns_motif_2), dtype=int)

    # Initalization maginales of T table
    for i in range(rows_motif_1):
        for j in range(columns_motif_1):
            for k in range(rows_motif_2):
                for l in range(columns_motif_2):
                    t_matrix[i][j][k][l] = 0
                    t_matrix[0][j][k][l] = (k + 1) * (l + 1)
                    t_matrix[i][0][k][l] = (k + 1) * (l + 1)
                    t_matrix[i][j][0][l] = (i + 1) * (j + 1)
                    t_matrix[i][j][k][0] = (i + 1) * (j + 1)

    # Filling of T table
    for i in range(1, rows_motif_1):
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

    treatment_time = time.time() - start_treatment
    return pre_treatment_time, treatment_time, t_matrix


rows = args.rows[0]
columns = args.columns[0]

motif_1 = utils.get_motif(rows, columns)
motif_2 = utils.get_motif(rows, columns)

# print(motif_1)
# print(motif_2)

start = time.time()
result = t_matrice_func(motif_1, motif_2)
end = time.time() - start

dtf = [result[0], result[1], end, rows, columns]
stats_file_name = "stats/sequential"

stat_file_name = 'stats/sequential/stats_for_seq_rows_{}_columns_{}.csv'.format(rows, columns)
utils.list_to_same_file(stat_file_name, [dtf])

print("Computations = ", result)
print("Execution time = ", end)
