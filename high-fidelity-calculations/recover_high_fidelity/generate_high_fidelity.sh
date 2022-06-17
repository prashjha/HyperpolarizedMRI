#!/bin/bash
MY_PWD=$(pwd)

OSTYPE="$(uname -s)"

if [[ $OSTYPE == "Darwin" ]]; then
	source ~/.bash_profile
	source ${HOME}/opt/anaconda3/etc/profile.d/conda.sh
else
	source ~/.bashrc
	source ${HOME}/anaconda3/etc/profile.d/conda.sh
fi

cd $MY_PWD

num_unc="3"
num_scans="30"

# oed dir
oed_dir="SNR_2_NGauss_5_NumberUncertain_3_Nscans_30_optim_FA_const_design_update_kpl_all_domain"
echo "oed_dir = $oed_dir"

tar_dir="${MY_PWD}/../oed_calculations/results/${oed_dir}/"
mkdir -p $tar_dir

echo "****** Directory = $tar_dir ******"
echo " "
cd $tar_dir
cp ${MY_PWD}/run_high_fidelity.py .
cp ${MY_PWD}/resamplefem.py .

# run high-fidelity model
if [[ $OSTYPE == "Darwin" ]]; then
	source ${HOME}/opt/anaconda3/etc/profile.d/conda.sh
else
	source ${HOME}/anaconda3/etc/profile.d/conda.sh
fi

# solve high-fidelity
conda activate confen
mpirun -n 24 python run_high_fidelity.py 

# resample to generate voxel data
if [[ $OSTYPE == "Darwin" ]]; then
	conda activate hpmri
	python resamplefem.py
fi
