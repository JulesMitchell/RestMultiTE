import pickle
from idtxl.stats import network_fdr

# Set target numbers to combine
subjects = ['sub-01']
timepoints = ['ses-01']
targets = [range(32)]

# Load results using pickle
for subject in subjects:
    for timepoint in timepoints:
        res_list = []
        for target_id in targets:
            path = 'my_directory/res.{}.{}.{}.pkl'.format(str(subject), str(timepoint), str(target_id))
            res_list.append(pickle.load(open(path, 'rb')))

        res = network_fdr({'alpha_fdr': 0.05}, *res_list)

        pickle.dump(res, open('result_{}.{}.p'.format(str(subject), str(timepoint)), 'wb'))