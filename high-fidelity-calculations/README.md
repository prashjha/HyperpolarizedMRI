# Scripts and files
## HPModel.m
Script provides a class implementing pharmacokinetic model for HP-MRI.

## def_params.mat
MATLAB data that hosts default model and design parameters and HP signals computed at default model and design parameters. 

# Directories

## calculate_reference_noise
In this folder, script is provided to computed the HP signals at default parameters. The script saves the data in `def_params.mat` that is described above.

## mesh
In this directory, one of the mesh of brain tissue along with vascular subdomain data is stored. 

## oed_calculations
In this directory, OED results are stored. 

In `oed_calculations/results/SNR_2_NGauss_5_NumberUncertain_3_Nscans_30_optim_FA_const_design_update_kpl_all_domain`, all results of high-fidelity simulation and inverse problems are stored. For high-fidelity simulation, we have only provided `nii` files and not the VTK files to minimize the memory. The `nii` files are required for inverse problems. 

The file `recover_all_snr_all_voxels.npz` inside `oed_calculations/results/SNR_2_NGauss_5_NumberUncertain_3_Nscans_30_optim_FA_const_design_update_kpl_all_domain/HF_recover_NumRecoverParams_1_NumTrials_25_Rel_Noise_1/`, stores recovered parameters for five SNR values (`SNR_data`) and for 25 selected cells. This file can be read in jupyter notebook `recover_high_fidelity/OED_plot.ipynb` to plot the results. 

## recover_high_fidelity
In this folder, bash and python scripts are provided to compute the high-fidelity model and solve inverse problems for parameter recovery. Further details can be found `README.md` file inside this folder.
