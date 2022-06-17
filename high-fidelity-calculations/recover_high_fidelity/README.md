# High-fidelity simulation

Script `environment.yml` can be used to create a python environment to run and test the codes. To create a conda environment usign this script, run following on terminal
```sh
conda env create -f environment.yml
```

## Scripts

1. `hp_model_class.py` - Implements high-fidelity HP-MRI model using Fenics

2. `resamplefem.py` -  Uses [c3d](http://www.itksnap.org/pmwiki/pmwiki.php?n=Convert3D.Documentation) to project finite element mesh on brain to solution on uniform mesh over a cube (cube containing original brain domain) and saves the data as `nii` format. The `nii` file provides HP signal for each cell of uniform mesh on cube.

3. `run_high_fidelity.py` - Python script that loads finite element mesh of a brain, creates finite element function space, instantiates HPModel class, and simulates the high-fidelity model. 

3. `generate_high_fidelity.sh` - Bash script that calls the python script `run_high_fidelity.py`. 

## Visualizing high-fidelity data and selecting cells for inverse problem
In jupyter notebook `plot_HF_data_on_voxels.ipynb`, high-fidelity data at different cells on grid are visualized. First 30 cells were identified based on **sign distance** of cells from vascular subdomain. For these 30 cells, the peak vascular pyruvate is used to sort them in decreasing peak vascular pyruvate. Out of 30, 25 cells are selected in four ranges of vascular pyruvate. The resulting sorted cells along with their cell (x,y,z) location are stored in numpy data `new_voxel_sorting_data.npz`.

# Solving inverse problem on high-fidelity data

For HP-signal data of selected 25 cells, we solve inverse problem to recover uncertain parameters. For noisy data, 5 different SNRs, `[2, 5, 10, 15, 20`, were considered. 

There are two scripts for this task: a MATLAB script that solves the optimization problem and a bash script that calls the MATLAB script.

Inverse problem can be tested as we have provided high-fidelity data in the directory `oed_calculations/results/SNR_2_NGauss_5_NumberUncertain_3_Nscans_30_optim_FA_const_design_update_kpl_all_domain/0/`. 

## Scripts

1. `recover_voxel.m` - Given a voxel id and SNR data among some other parameters, this scripts first reads HP-signal data from high-fidelity simulation and creates samples of noisy data. Next, it takes `NumTrials` number of samples of noisy data and solves the inverse problem to determine the uncertain parameters. 

2. `run_voxel.sh` - Bash script that calls `recover_voxel.m` based on provided input values. Since we have 25 cells to consider for inverse problem, we divide these 25 cells into 5 chunks, e.g., chunk 1 has `[1,2,3,4,5]`, chunk 2 has `[,6,7,8,9,10]`, and so on. From input, we can select one of the chunk and then this script will run inverse problem in parallel for 5 cells in that given chunk. 

# Plotting OED results
Jupyter notebook `OED_plot.ipynb` reads and plots the OED results for different cells (25 cells were considered) for a given `SNR_data` value.

