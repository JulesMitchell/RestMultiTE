# te_hpc.py 

"""
The following code file was called by a PBS script (run_network.pbs) passing the following:
1. Target index (i.e. target channel)
2. File path to csv data files. 
3. Subject ID
4. Session number
5. Resting state condition

After loading, the channel column is removed and the time-series is passed to the IDTxl data object with dimension order of processes (i.e. channels) and samples. 
"""
import sys
import os
from pathlib import Path
import pickle
import numpy as np
import pandas as pd

sys.path.append('/new_idtxl')

from idtxl.multivariate_te import MultivariateTE
from idtxl.data import Data
from idtxl import idtxl_io as io
import idtxl.estimators_jidt

OUTPUT_PATH = Path('/data/pickles')

if not OUTPUT_PATH.exists():
    print(f"Your output path doesnt exist {OUTPUT_PATH}")

# Read parameters from shell call
target_id = int(sys.argv[1]) # Target index
file_path = str(sys.argv[2]) # File directory
subject_id = sys.argv[3] # Subject 
session_id = sys.argv[4] # Session
condition = sys.argv[5]

# Load time series from csv
time_series = pd.read_csv(file_path, header=None) # adjust as needed
time_series = time_series.drop(columns = 0)
#time_series = time_series.iloc[:, 10000:15000] # Cut the data to speed up analysis during testing

# Initialise Data object and set dim_order to reflect your data
dat = Data(time_series, dim_order='ps')

# Initialise analysis object and define settings
network_analysis = MultivariateTE()
settings = {'cmi_estimator': 'JidtGaussianCMI', #OpenCLKraskovCMI/JidtGaussianCMI/JidtKraskovCMI
            'max_lag_sources': 5,
            'min_lag_sources': 1,
            'n_perm_max_stat': 200,
            'alpha_max_stat': 0.05,
            'n_perm_min_stat': 200,
            'alpha_min_stat': 0.05,
            'permute_in_time': True,
            'perm_type': 'random',
            'n_perm_omnibus': 650,
            'alpha_omnibus': 0.05,
            'n_perm_max_seq': 650,
            'alpha_max_seq': 0.05
}

# Run analysis
res = network_analysis.analyse_single_target(settings, dat, target_id)

# Construct the output file path
file_path = OUTPUT_PATH / f'{subject_id}'/f'{session_id}'/f'{subject_id}_{session_id}_{condition}_#{target_id}.p'

# File save
with open(file_path, 'wb') as f:
    pickle.dump(res, f)
