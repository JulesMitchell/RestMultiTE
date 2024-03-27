# Conda environment installation
# git clone https://github.com/pwollstadt/IDTxl.git
conda create --name idtxl python=3.11 pip matplotlib h5py scipy networkx
conda activate idtxl
conda install -c conda-forge jpype1  # required by CPU JIDT estimators
conda install -c conda-forge pyopencl  # required by GPU OpenCL estimators

cd IDTxl
pip install -e .
python demos/demo_bivariate_mi.py