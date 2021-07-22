import utils
import serial
import time


# This is a sample Python script.

# Press Shift+F10 to execute it or replace it with your code.
# Press Double Shift to search everywhere for classes, files, tool windows, actions, and settings.


def print_hi(name):
    # Use a breakpoint in the code line below to debug your script.
    print(f'Hi, {name}')  # Press Ctrl+F8 to toggle the breakpoint.


# Press the green button in the gutter to run the script.
if __name__ == '__main__':
    rows = 50
    columns = 50

    motif_1 = utils.get_motif(5, columns)
    motif_2 = utils.get_motif(rows, columns)

    print(motif_1)
    print(motif_2)

    start = time.time()

    print("Pre-treatment computing time = ", serial.t_matrice_func(motif_1, motif_2))

    end = time.time() - start
    print("Computation time = ", end)

# See PyCharm help at https://www.jetbrains.com/help/pycharm/
