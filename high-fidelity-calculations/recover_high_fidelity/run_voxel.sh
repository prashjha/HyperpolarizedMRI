#!/bin/bash
MY_PWD=$(pwd)

OSTYPE="$(uname -s)"
WHICHPC="HomeDesktop"
MATLABEXE="${HOME}/.local/MATLAB/R2020b/bin/matlab"
if [[ $OSTYPE == "Darwin" ]]; then
  source ~/.bash_profile
  WHICHPC="Macbook"
  MATLABEXE="/Applications/MATLAB_R2020b.app/bin/matlab"
else
  source ~/.bashrc
  if [[ -f "${HOME}/UTDesktop.txt" ]]; then
		WHICHPC="UTDesktop"
		MATLABEXE="${HOME}/.local/MATLAB/R2021a/bin/matlab"
	else
		WHICHPC="Desktop"
		MATLABEXE="${HOME}/.local/MATLAB/R2022a/bin/matlab"
	fi
fi

echo "OSTYPE = $OSTYPE, WHICHPC = $WHICHPC, MATLABEXE = $MATLABEXE"

# SNR data to consider
declare -a SNR_data_vec=("2" "5" "10" "15" "20")

## collect input arguments
# input 1: SNR to generate noisy data
# $1 \in [0, 1, 2, 3, 4]
# input 2: Number of recovered parameters
# $2 \in [1, 2, 3]
# input 3: Recover true (0 - false, 1 - true)
# $3 \in [0, 1]
# input 4: Voxel chunk to consider
# $4 \in [0, 1, 2, 3, 4]
SNR_data="${SNR_data_vec[$1]}"
num_rec_param="$2"
rec_true="$3"

voxel_chunk="$4"
voxel_lw="$(expr ${voxel_chunk} \* 5 + 1)"
voxel_up="$(expr ${voxel_lw} + 4)"

echo "SNR_data = ${SNR_data}, voxel_chunk = ${voxel_chunk}, voxel_lw = ${voxel_lw}, voxel_up = ${voxel_up}"

num_trials="25"
num_scans="30"
rel_noise="1"

# loop over voxels (divide 25 voxels in chunks of 5)
for voxel_i in $(seq ${voxel_lw} 1 ${voxel_up}); do

	cd $MY_PWD
	echo "voxel_i = ${voxel_i}"

	# oed dir
	oed_dir="SNR_2_NGauss_5_NumberUncertain_3_Nscans_30_optim_FA_const_design_update_kpl_all_domain"
	echo "oed_dir = $oed_dir"

	# create directory within oed dir to store inverse problem results
	rec_dir=`printf "HF_recover_NumRecoverParams_%d_NumTrials_%d_Rel_Noise_%d/SNR_data_%d/voxel_%d" ${num_rec_param} ${num_trials} ${rel_noise} ${SNR_data} ${voxel_i}`
	if [[ ${rec_true} == "1" ]]; then
		rec_dir=`printf "HF_recover_NumRecoverParams_%d_Rec_True/No_noise/voxel_%d" ${num_rec_param} ${voxel_i}`
	fi	
	echo "rec_dir = $rec_dir"

	tar_dir="${MY_PWD}/../oed_calculations/results/${oed_dir}/${rec_dir}/"
	mkdir -p $tar_dir
	cp recover_voxel.m $tar_dir

	cd $tar_dir
	echo "pwd = $(pwd)"

	# screen name
	scr_nm=`printf "HPMRI-%d-%d" ${SNR_data} ${voxel_i}`
	echo "scr_nm = ${scr_nm}"

	# run
	if ! screen -list | grep -q "${scr_nm}"; then
		echo "Running script"
    screen 	-dmS "${scr_nm}" \
					${MATLABEXE} -nodesktop -r \
					"a = recover_voxel(${num_rec_param}, ${num_trials}, ${SNR_data}, ${rec_true}, ${rel_noise}, ${voxel_i}); exit()"
	fi
	
	# ${MATLABEXE} -nodesktop -r \
	# 				"a = recover_voxel(${num_rec_param}, ${num_trials}, ${SNR_data}, ${voxel_i}); exit()"

done
