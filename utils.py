import csv

ALL_DATASET_PATH = "input/datasets.txt"


def clean_str(label):
    return label[:len(label) - 1]


def get_dataset_path(dataset_pos=0):
    with open(ALL_DATASET_PATH, 'r') as all_dataset_paths_file:
        lines = all_dataset_paths_file.readlines()

    return clean_str(str(lines[dataset_pos]))


def get_dataset(file_path):
    with open(file_path, 'r') as file:
        dataset = file.readline()

    return dataset


def get_motif(rows=10, columns=10, dataset_pos=0):
    dataset_path = get_dataset_path(dataset_pos)
    dataset = get_dataset(dataset_path)
    # print(dataset[:50])
    motif = []
    step = 0
    for cpt in range(rows):
        motif.append(dataset[step:step + columns])
        step += columns

    return motif


# 1D Sequential alignment function
def classic_seq_edit_distance(motif_1, motif_2):
    assert isinstance(motif_1, str)
    assert isinstance(motif_2, str)

    motif_1_length = len(motif_1)
    motif_2_length = len(motif_2)

    scores = []
    row = motif_1_length + 1
    column = motif_2_length + 1

    # Initializing the scores list
    for i in range(row):
        inter = []
        for j in range(column):
            inter.append(0)
        scores.append(inter)

    for i in range(row):
        scores[i][0] = i

    for j in range(column):
        scores[0][j] = j

    for i in range(1, row):
        for j in range(1, column):
            if motif_1[i - 1] == motif_2[j - 1]:
                scores[i][j] = scores[i - 1][j - 1]
                continue
            x = scores[i - 1][j] + 1  # Deletion
            y = scores[i][j - 1] + 1  # Insertion
            z = scores[i - 1][j - 1] + 1  # Substitution

            scores[i][j] = min(x, y, z)

    result = scores[motif_1_length][motif_2_length]
    # del scores

    return scores


def del_func(char):
    return 1


def ins_func(char):
    return 1


# filling of the Dr Table
def dr_matrice_func(x):
    motif = x
    assert isinstance(motif, list)
    assert isinstance(motif[0], str)

    dr_matrix = []

    row = len(motif)
    column = len(motif[0])

    # Initializing the scores list
    for i in range(row):
        inter = []
        for j in range(column):
            inter.append(0)
        dr_matrix.append(inter)

    for i in range(0, row):
        for j in range(0, column):
            step = 0
            for p in range(j):
                substring = motif[i][p]
                step += del_func(substring)
            dr_matrix[i][j] = step

    return dr_matrix


# filling of the Dc Table
def dc_matrice_func(x):
    motif = x
    assert isinstance(motif, list)
    assert isinstance(motif[0], str)

    dc_matrix = []

    row = len(motif)
    column = len(motif[0])

    # Initializing the scores list
    for i in range(row):
        inter = []
        for j in range(column):
            inter.append(0)
        dc_matrix.append(inter)

    for i in range(0, row):
        for j in range(0, column):
            step = 0
            for p in range(i):
                substring = motif[p][j]
                step += del_func(substring)
            dc_matrix[i][j] = step

    return dc_matrix


# filling of the Ir Table
def ir_matrice_func(y):
    motif = y
    assert isinstance(motif, list)
    assert isinstance(motif[0], str)

    ir_matrix = []

    row = len(motif)
    column = len(motif[0])

    # Initializing the scores list
    for i in range(row):
        inter = []
        for j in range(column):
            inter.append(0)

        ir_matrix.append(inter)

    for i in range(0, row):
        for j in range(0, column):
            step = 0
            for p in range(j):
                substring = motif[i][p]
                step += ins_func(substring)
            ir_matrix[i][j] = step

    return ir_matrix


# filling of the Ic Table
def ic_matrice_func(y):
    motif = y
    assert isinstance(motif, list)
    assert isinstance(motif[0], str)

    ic_matrix = []

    row = len(motif)
    column = len(motif[0])

    # Initializing the scores list
    for i in range(row):
        inter = []
        for j in range(column):
            inter.append(0)
        ic_matrix.append(inter)

    for i in range(0, row):
        for j in range(0, column):
            step = 0
            for p in range(i):
                substring = motif[p][j]
                step += del_func(substring)
            ic_matrix[i][j] = step

    return ic_matrix


# Filling of Dc table
# def dc_matrice_func(rows=10, columns=10):
#     dc_matrix = []
#     for i in range(rows):
#         inter = []
#         for j in range(columns):
#             inter.append(0)
#         dc_matrix.append(inter)
#
#     for i in range(rows):
#         for j in range(columns):
#             dc_matrix[i][j] = i + 1
#
#     return dc_matrix


# Filling of R table
def r_matrice_func(x, y):
    r_matrix = []

    assert isinstance(x, list)
    assert isinstance(y, list)

    dr = dr_matrice_func(x)
    ir = ir_matrice_func(y)

    row_1 = len(x)
    col_1 = len(x[0])
    row_2 = len(y)
    col_2 = len(y[0])

    # 4 dimension matrix initializations
    for i in range(row_1):
        inter_1 = []
        for j in range(col_1):
            inter_2 = []
            for k in range(row_2):
                inter_3 = []
                for l in range(col_2):
                    inter_3.append(0)
                inter_2.append(inter_3)
            inter_1.append(inter_2)
        r_matrix.append(inter_1)

    # Building result
    for i in range(row_1):
        for j in range(col_1):
            for k in range(row_2):
                for l in range(col_1):
                    r_matrix[i][j][k][l] = dr[i][j] + ir[k][l]
                    # r_matrice[i][j][k][l] = classic_seq_edit_distance(motif_1[i][0:j], motif_2[k][0:l])

    return r_matrix


# Filling of C table
def c_matrice_func(x, y):
    c_matrix = []

    assert isinstance(x, list)
    assert isinstance(y, list)

    dc = dc_matrice_func(x)
    ic = ic_matrice_func(y)

    row_1 = len(x)
    col_1 = len(x[0])
    row_2 = len(y)
    col_2 = len(y[0])

    # 4 dimension matrix initializations
    for i in range(row_1):
        inter_1 = []
        for j in range(col_1):
            inter_2 = []
            for k in range(row_2):
                inter_3 = []
                for l in range(col_2):
                    inter_3.append(0)
                inter_2.append(inter_3)
            inter_1.append(inter_2)
        c_matrix.append(inter_1)

    # Building result
    for i in range(row_1):
        for j in range(col_1):
            for k in range(row_2):
                for l in range(col_1):
                    c_matrix[i][j][k][l] = dc[i][j] + ic[k][l]
                    # r_matrice[i][j][k][l] = classic_seq_edit_distance(motif_1[i][0:j], motif_2[k][0:l])

    return c_matrix


def split_data(data, nprocs, rank):
    r, diff = divmod(len(data), nprocs)
    counts = [r + 1 if p < diff else r for p in range(nprocs)]

    # determine the starting and ending indices of each sub-task
    starts = [sum(counts[:p]) for p in range(nprocs)]
    ends = [sum(counts[:p + 1]) for p in range(nprocs)]

    # converts data into a list of arrays
    datas = [data[starts[p]:ends[p]] for p in range(nprocs)]

    return datas[rank]


def list_to_same_file(file, datas):
    with open(file, 'a') as k_anon_table:
        writer = csv.writer(k_anon_table, delimiter=';')
        for item in datas:
            writer.writerow(item)

# print(classic_seq_edit_distance("vadele", "vadelo"))

# print(get_motif(5, 5))
# print(get_motif(5, 5, 1))

# print(r_matrice_func(get_motif(100, 100), get_motif(100, 100, 1))[:1][0][0][0][0])
# print(c_matrice_func(get_motif(100, 100), get_motif(100, 100, 1))[:1][0][0][0][0])
# print(t_matrice_func(get_motif(50, 50), get_motif(50, 50, 1))[:1][0][0][0][0])
