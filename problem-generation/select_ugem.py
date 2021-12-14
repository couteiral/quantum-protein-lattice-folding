"""
This program is used to select the
UGEM (unique ground energy minimum)
structures from the results of the
calculations.

Carlos Outeiral
University of Oxford
June 2019
"""

import os
import numpy as np
from tqdm import tqdm
from itertools import product


with open('mj_models', 'r') as mj_file:
    seqs = mj_file.readlines()
    seqs = [seq[:-1] for seq in seqs]

for length in tqdm(range(6, 10), desc='Screening all MJ models...'):

    # Gather all histograms
    seq_list = [seq for seq in seqs if len(seq) == length]
    hist_list = []
    prefix = 'lat_2D_s%d' % length
    for seq in tqdm(seq_list, 
                    desc='Checking 2D models of size %d' % length,
                    total=len(seq_list)):

        new_prefix = seq + '_' + prefix
        hist = np.loadtxt(new_prefix + '.dist').reshape(-1, 2)
        hist_list.append(hist)

    # Find what is the minimum possible number of structures
    cur_min = 8
    for hist in hist_list:
        if hist[-1, 1] < cur_min:
            cur_min = hist[-1, 1]

    ugem_list = []
    # Gather structures with this number of structures
    for i, hist in enumerate(hist_list):
        if hist.shape[0] == 1:
            continue
        if hist[-1, 1] == cur_min:
            ugem_list.append(seq_list[i] + '\n')

    print(ugem_list)

    with open(prefix + '_MJ.ugem', 'w') as f:
        f.writelines(ugem_list)

