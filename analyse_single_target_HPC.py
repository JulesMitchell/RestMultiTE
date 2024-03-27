# analyse_single_target.py

import sys
import os
import pickle
import numpy as np
import pandas as pd

from idtxl.multivariate_te import MultivariateTE
from idtxl.data import Data
from idtxl import idtxl_io as io
import idtxl.estimators_jidt

# Read parameters from shell call
target_id = int(sys.argv[1]) # Target index
file_dir = sys.argv[2] # File directory
subject_id = sys.argv[3] # Subject 
session_id = sys.argv[4] # Session

# Change to the subject-session file directory
os.chdir(file_dir)

# Load time series from numpy array
time_series = pd.load('EO_chanactivity_.csv') # adjust as needed
time_series = time_series.drop(columns = 0)

# Initialise Data object and set dim_order to reflect your data
dat = Data(time_series, dim_order='psr')

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
            'n_perm_omnibus': 500,
            'alpha_omnibus': 0.05,
            'n_perm_max_seq': 500,
            'alpha_max_seq': 0.05
}

# Run analysis
res = network_analysis.analyse_single_target(settings, dat, target_id)

# Construct the output file path
output_dir = os.path.join('derivatives', subject_id, session_id) # add project specific folder
os.makedirs(output_dir, exist_ok=True)  # Create the directory if it doesn't exist
file_path = os.path.join(output_dir, f'{subject_id}_{session_id}_{target_id}_EO_target.p')

# File save
with open(file_path, 'wb') as f:
    pickle.dump(res, f)