import utils
import time
import argparse
import sys

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
    t_matrice = []
    # 4 dimension matrix initializations
    for i in range(rows_motif_1):
        assert isinstance(motif_1[0], str)
        inter_1 = []
        for j in range(columns_motif_1):
            inter_2 = []
            for k in range(rows_motif_2):
                assert isinstance(motif_2[0], str)
                inter_3 = []
                for l in range(columns_motif_2):
                    inter_3.append(0)
                inter_2.append(inter_3)
            inter_1.append(inter_2)
        t_matrice.append(inter_1)

    # Initalization maginales of T table
    for i in range(rows_motif_1):
        for j in range(columns_motif_1):
            for k in range(rows_motif_2):
                for l in range(columns_motif_2):
                    t_matrice[i][j][k][l] = 0
                    t_matrice[0][j][k][l] = (k + 1) * (l + 1)
                    t_matrice[i][0][k][l] = (k + 1) * (l + 1)
                    t_matrice[i][j][0][l] = (i + 1) * (j + 1)
                    t_matrice[i][j][k][0] = (i + 1) * (j + 1)

    # Filling of T table
    for i in range(1, rows_motif_1):
        for j in range(1, columns_motif_1):
            for k in range(1, rows_motif_2):
                for l in range(1, columns_motif_2):
                    val1 = max(t_matrice[i - 1][j][k][l] + dr_mat[i][rows_motif_1 - 1],
                               t_matrice[i][j - 1][k][l] + dc_mat[i][rows_motif_1 - 1])
                    val2 = max(val1, t_matrice[i][j][k - 1][l] + ir_mat[i][rows_motif_2 - 1])
                    val3 = max(val2, t_matrice[i][j][k - 1][l - 1] + ic_mat[i][rows_motif_2 - 1])
                    val4 = max(val3, t_matrice[i - 1][j][k - 1][l] + r_mat[i][j][k][l])
                    val5 = max(val4, t_matrice[i][j - 1][k][l - 1] + c_mat[i][j][k][l])
                    val6 = max(val5, t_matrice[i - 1][j - 1][k - 1][l - 1] + c_mat[i - 1][j][k - 1][l] +
                               r_mat[i][j][k][l])
                    t_matrice[i][j][k - 1][l] = max(val6, t_matrice[i - 1][j - 1][k - 1][l - 1] + c_mat[i][j][k][l] +
                                                    r_mat[i][j - 1][k][l - 1])

    # for i in range(rows_motif_1):
    #     for j in range(columns_motif_1):
    #         for k in range(rows_motif_2):
    #             for l in range(columns_motif_2):
    #                 print("[{}][{}][{}][{}] = ".format(i, j, k, l), t_matrice[i][j][k][l])

    treatment_time = time.time() - start_treatment
    return pre_treatment_time, treatment_time


rows = args.rows[0]
columns = args.columns[0]

motif_1 = utils.get_motif(rows, columns)
motif_2 = utils.get_motif(rows, columns)

print(motif_1)
print(motif_2)

start = time.time()
result = t_matrice_func(motif_1, motif_2)
end = time.time() - start

dtf = [result[0], result[1], end, rows, columns]
stats_file_name = "stats/sequential"

stat_file_name = 'stats/sequential/stats_for_seq_rows_{}_columns_{}.csv'.format(rows, columns)
utils.list_to_same_file(stat_file_name, [dtf])

print("Computations = ", result)
print("Execution time = ", end)
