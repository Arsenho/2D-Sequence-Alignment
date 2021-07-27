import os

PROC_RANGE = [2, 4, 8, 16, 32, 50]
rows = 120
columns = 120

for proc in PROC_RANGE:
    cmd = "mpirun -n {} python3 parallel.py --rows {} --columns {}".format(proc, rows, columns)
    os.system(cmd)
