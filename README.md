Code (R, Python) and BRMS supplementary output for analysis of EEG resting state multivariate transfer entropy and network topology prior to and following oral low-dose ketamine treatment. Note, the EEG files were pre-processed in EEGLAB - code found here: https://github.com/JulesMitchell/RestComplexity

Order of code applied:

Python: Processed data (.csv files) analaysed with IDTxL package using 'te_entropy.py'. Data (directed acyclic graph) exported as pickle file. Network metrics extracted with 'network_feature_extraction.py'.
R: Data analysed using BRMS (STAN backend) mixed effect models. Code files structured by metric. 
Python: Posterior effects plotted using 'topoplots.py'.

You are free to use this or any other code from this repository for your own projects and publications. Citation or reference to the repository is not required, but would be much appreciated (see more on README.md).

Any code feedback is welcome.

Code created by Jules Mitchell March 2025.
