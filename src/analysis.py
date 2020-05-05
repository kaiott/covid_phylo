import align_tools as at
import matplotlib.pyplot as plt
import numpy as np
from collections import Counter

def h(x):
    if x>0:
        return 1
    else:
        return 0

def get_counter(arr, lower_sat=None, upper_sat=None):
    result = {}
    for val in arr:
        if (upper_sat is None or val < upper_sat) and (lower_sat is None or val > lower_sat):
            result[val] = result.get(val, 0) + 1
        elif upper_sat is not None and val >= upper_sat:
            result[upper_sat] = result.get(upper_sat, 0) + 1
        else:
            result[lower_sat] = result.get(lower_sat, 0) + 1

    return result


def analyse_gaps(num_gaps, collaps_factor=1):
    print(get_counter(num_gaps, 1))
    has_gaps = [h(num_gap) for num_gap in num_gaps]
    num_gaps_collaps = [sum(has_gaps[max([collaps_factor*i, 0]):min([collaps_factor*(i+1), len(has_gaps)])]) for i in range(int(len(has_gaps)/collaps_factor)+1)]

    ax = plt.subplot(111)
    x = [n for n in num_gaps_collaps]
    ax.bar(range(len(num_gaps_collaps)), num_gaps_collaps)
    plt.xlabel('x')
    plt.ylabel('y')
    plt.title('Posiciones con gaps')
    plt.show()


def analyse_changes(num_vars_det, num_vars_all):
    vars_det_sites = get_counter(num_vars_det, 0, 4)
    vars_all_sites = get_counter(num_vars_all, 0, 4)
    print('only determined')
    print([f'k={k}: {vars_det_sites.get(k, 0)}, {vars_det_sites.get(k, 0) / len(num_vars_det) * 100:.2f}%' for k in vars_det_sites])
    print('also undetermined')
    print([f'k={k}: {vars_all_sites.get(k, 0)}, {vars_all_sites.get(k, 0) / len(num_vars_all) * 100:.2f}%' for k in vars_all_sites])

    x = [n for n in vars_det_sites]
    y = [vars_det_sites.get(n, 0) for n in x]
    z = [vars_all_sites[n] for n in x]

    ax = plt.subplot(111)
    bar1 = ax.bar(np.array(x)-0.1, y, width=0.2, color='b', align='center')
    bar2 = ax.bar(np.array(x)+0.1, z, width=0.2, color='r', align='center')
    ax.legend( (bar1[0], bar2[0]), ('Solo bases conocidas', 'Incluyendo bases desconocidas'))
    plt.xlabel('k (saturación en 4)')
    plt.xticks([1, 2, 3, 4])
    plt.ylabel('n_k')
    plt.title('Histograma de nucleotidos distintos por posición')
    plt.show()



records = at.aligned_records_by_tag("complete")
num_gaps, num_vars_det, num_vars_all = at.analyse_alignment(records)
print("done anaylsis")

analyse_gaps(num_gaps, collaps_factor=300)

analyse_changes(num_vars_det, num_vars_all)
