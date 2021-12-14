"""

This script is used to couple the
programs contact_bf and energy_calc,
to obtain a list of really hard
lattice protein problems.

Carlos Outeiral
University of Oxford
May 2019
"""

import os
import numpy as np
from tqdm import tqdm
from itertools import product

np.random.seed(448)
mj_model = ['Y', 'C', 'K', 'E', 'L', 'V', 'R', 'H',
            'P', 'M', 'A', 'W', 'S', 'I', 'T', 'Q',
            'F', 'D', 'G', 'N']


# First, we build all the contact lists
for dim in [2, 3]:

    if dim == 2:
        for length in tqdm(range(4, 16), desc='Building %dD lists of contacts...' % dim, total=12):
            os.system('./contact_bf %d %d lat_%dD_s%d >> log.log' % (length, dim, dim, length))
    if dim == 3:
        for length in tqdm(range(4, 12), desc='Building %dD lists of contacts...' % dim, total=8):
            os.system('./contact_bf %d %d lat_%dD_s%d >> log.log' % (length, dim, dim, length))

with open('run_script.sh', 'w') as run_script:
    for length in tqdm(range(6, 10), desc='Sampling 2D MJ models...'):

        drawn_samples = []
        prefix = 'lat_2D_s%d' % length
        n_samples = 5000

        for _ in tqdm(range(n_samples), desc='Testing models of size %d' % length):

            # Two-dimensional
            seq = ''.join(np.random.choice(mj_model, size=length, replace=True))
            new_prefix = seq + '_' + prefix
            if seq in drawn_samples or seq[::-1] in drawn_samples:
                continue
            else:
                drawn_samples.append(seq)

            with open(seq + '_' + prefix + '.sh', 'w') as afile:
                afile.write('./energy_calc %s %s MJ > /dev/null\n' % (prefix, seq))
                afile.write('python analyse_sols.py %s\n' % new_prefix)
                afile.write('rm %s.sol\n' % new_prefix)

            run_script.write('bash %s.sh\n' % new_prefix)

    for length in tqdm(range(6, 10), desc='Sampling 3D MJ models...'):

        drawn_samples = []
        prefix = 'lat_3D_s%d' % length
        n_samples = 5000

        for _ in tqdm(range(n_samples), desc='Testing models of size %d' % length):

            seq = ''.join(np.random.choice(mj_model, size=length, replace=True))
            new_prefix = seq + '_' + prefix
            if seq in drawn_samples or seq[::-1] in drawn_samples:
                continue
            else:
                drawn_samples.append(seq)

            with open(seq + '_' + prefix + '.sh', 'w') as afile:
                afile.write('./energy_calc %s %s MJ > /dev/null\n' % (prefix, seq))
                afile.write('python analyse_sols.py %s\n' % new_prefix)
                afile.write('rm %s.sol\n' % new_prefix)

            run_script.write('bash %s.sh\n' % new_prefix)
