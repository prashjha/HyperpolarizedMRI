function rec_params = recover_voxel(NumberRecParams, NumTrials, SNRData, RecTrue, RelativeNoise, Voxel_ID)

disp('recover_voxel')
NumberRecParams, NumTrials, SNRData, Voxel_ID, RelativeNoise

% generate different sequences for different SNR data
if SNRData < 0
	rng(0)
else
	rng(SNRData)
end

% Nscans
Nscans = 30;

% find reference peak Mxy
def_params = load('../../../../../../def_params.mat')
def_params = def_params.str;
Mxy_pyr_max = def_params.Mxy_pyr_max;

% noise in individual signal
Mxy_std = Mxy_pyr_max / SNRData;

% directories
OED_res_dir = '../../../';
rec_dir = './';
HF_sim_dir = '../../../0/';

fname_log = sprintf('%s%s', rec_dir, '/recover.log');
fname_data = sprintf('%s%s', rec_dir, '/recover.mat');
diary(fname_log);

% add paths
addpath ../../../.
addpath ../../../../../../.

% load synthetic data
mymask = niftiread(sprintf('%spyrlacmip.nii', HF_sim_dir));
Nvoxel = 16
oedpyr = zeros(Nscans, Nvoxel, Nvoxel, Nvoxel);
oedlac = zeros(Nscans, Nvoxel, Nvoxel, Nvoxel);
oedpyr_vasc = zeros(Nscans, Nvoxel, Nvoxel, Nvoxel);
for iddata = 1:Nscans
   disp(sprintf('iddata = %d',iddata));
   f = sprintf('%spyruvate%06d.nii', HF_sim_dir, iddata-1)
   oedpyr(iddata,:,:,:) = niftiread(f);
   f = sprintf('%slactate%06d.nii', HF_sim_dir, iddata-1)
   oedlac(iddata,:,:,:) = niftiread(f);
   f = sprintf('%spyruvatevasc%06d.nii', HF_sim_dir, iddata-1)
   oedpyr_vasc(iddata,:,:,:) = niftiread(f);
end

% select voxels
xroi =  [9, 7, 15, 15, 12, 7, 8, 2, 4, 14, 13, 13, 9, 14, 9, 7, 13, 6, 8, 3, 5, 9, 10, 2, 8]';
yroi =  [4, 15, 8, 8, 12, 9, 13, 9, 7, 13, 5, 8, 10, 6, 9, 9, 14, 8, 2, 9, 5, 9, 6, 8, 13]';
zroi =  [14, 7, 6, 6, 11, 7, 7, 8, 9, 10, 8, 7, 12, 9, 12, 9, 6, 8, 11, 12, 8, 8, 7, 7, 10]';

% OED parameters
% use const design OED for 'model SNR = 2' --> FaP = 35, FaL = 28
OED_FaP = 35; OED_FaL = 28;

% create a model and generate pure ground thruth
model = HPModel();
M0 = [0; 0];

% model parameters same as default parameters!!
params = def_params.params; 

% modify OED parameters
params.FaList(1,:) = OED_FaP*pi/180.; 
params.FaList(2,:) = OED_FaL*pi/180.;

% get model parameters to save as reference
T1P_ref = def_params.params.T1s(1); T1L_ref = def_params.params.T1s(2);
kpl_ref = def_params.params.ExchangeTerms(1, 2);
klp_ref = def_params.params.ExchangeTerms(2, 1);
kve_ref = def_params.params.kve(1);
ve_ref = def_params.params.ve(1); % ve in params is a column vector of two elements
kve_ve_ref = kve_ref / ve_ref;
t0_ref = def_params.params.t0(1);
VIF_scale_ref = def_params.params.scaleFactor(1); % scaling of VIF source
gammaPdfA = def_params.params.gammaPdfA(1); gammaPdfB = def_params.params.gammaPdfB(1);

% compute Mxy from data (use OED flip angles)
Nroivoxel = length(xroi);
Mxy_data_voxels = zeros(2, Nscans, Nroivoxel, NumTrials+1);
Mxy_std_vec = zeros(2, Nroivoxel);
for jjj =1:Nroivoxel 
	for kkk = 1:Nscans
		k_pyr = oedpyr(kkk, xroi(jjj), yroi(jjj), zroi(jjj));
		k_lac = oedlac(kkk, xroi(jjj), yroi(jjj), zroi(jjj));
		k_pyr_v = oedpyr_vasc(kkk, xroi(jjj), yroi(jjj), zroi(jjj));
		k_lac_v = 0.;
		k_sin_p = sin(params.FaList(1, kkk));
		k_sin_l = sin(params.FaList(2, kkk));
		% Mxy = sin(theta) * (pyr + pyr_vasc)
  	Mxy_data_voxels(1, kkk, jjj, 1) = k_sin_p * (k_pyr + k_pyr_v);
  	Mxy_data_voxels(2, kkk, jjj, 1) = k_sin_l * (k_lac + k_lac_v);
  end
  % noise
	if RelativeNoise == 1
		Mxy_std_vec(1, jjj) = max(Mxy_data_voxels(1, :, jjj, 1)) / SNRData;
		Mxy_std_vec(2, jjj) = max(Mxy_data_voxels(1, :, jjj, 1)) / SNRData;
	else
		Mxy_std_vec(1, jjj) = Mxy_std;
		Mxy_std_vec(2, jjj) = Mxy_std;
	end
end

% create arrays to store recovered parameters
rec_params = zeros(NumTrials+1, NumberRecParams, Nroivoxel);

% xstar -- ground truth, x0 -- initial guess

% recovery (optimizing) variables
switch (NumberRecParams)
	case(1) 
		kpl = optimvar('kpl','LowerBound',0);
	  % ground truth 
	  xstar.kpl = kpl_ref;
	  % random initial guess
	  x0.kpl        = unifrnd(.06,.24);
	  % remaining params
	  T1P = T1P_ref; T1L = T1L_ref; kve_ve = kve_ve_ref; 
	  t0 = t0_ref; VIF_scale = VIF_scale_ref;
	case(2) 
		kpl = optimvar('kpl','LowerBound',0);
	  kve_ve = optimvar('kve_ve','LowerBound',0); 
	  % ground truth 
	  xstar.kpl = kpl_ref;
	  xstar.kve_ve = kve_ve_ref;
	  % random initial guess
	  x0.kpl        = unifrnd(.06,.24);
	  x0.kve_ve      = unifrnd(.02,10);
	  % remaining params
	  T1P = T1P_ref; T1L = T1L_ref; 
	  t0 = t0_ref; VIF_scale = VIF_scale_ref;
	case(3) 
		kpl = optimvar('kpl','LowerBound',0);
	  kve_ve = optimvar('kve_ve','LowerBound',0); 
	  t0 = optimvar('t0','LowerBound',0); 
	  % ground truth 
	  xstar.kpl = kpl_ref;
	  xstar.kve_ve = kve_ve_ref;
	  xstar.t0 = t0_ref;
	  % random initial guess
	  x0.kpl        = unifrnd(.06,.24);
	  x0.kve_ve      = unifrnd(.02,10);
	  x0.t0         = unifrnd(0  ,8);
	  % remaining params
	  T1P = T1P_ref; T1L = T1L_ref; 
	  VIF_scale = VIF_scale_ref;
	case(4) 
	  kpl = optimvar('kpl','LowerBound',0);
	  kve_ve = optimvar('kve_ve','LowerBound',0); 
	  T1P = optimvar('T1P','LowerBound',0); 
	  T1L = optimvar('T1L','LowerBound',0); 
	  % ground truth 
	  xstar.kpl = kpl_ref;
	  xstar.kve_ve = kve_ve_ref;
	  xstar.T1P = T1P_ref;
	  xstar.T1L = T1L_ref;
	  % random initial guess
	  x0.kpl        = unifrnd(.06,.24);
	  x0.kve_ve      = unifrnd(.02,10);
	  x0.T1P        = unifrnd(20 ,40);
	  x0.T1L        = unifrnd(15 ,35);
	  % remaining params
	  t0 = t0_ref; VIF_scale = VIF_scale_ref;
	case(5) 
	  kpl = optimvar('kpl','LowerBound',0);
	  kve_ve = optimvar('kve_ve','LowerBound',0); 
	  T1P = optimvar('T1P','LowerBound',0); 
	  T1L = optimvar('T1L','LowerBound',0); 
	  t0 = optimvar('t0','LowerBound',0); 
	  % ground truth 
	  xstar.kpl = kpl_ref;
	  xstar.kve_ve = kve_ve_ref;
	  xstar.T1P = T1P_ref;
	  xstar.T1L = T1L_ref;
	  xstar.t0 = t0_ref;
	  % random initial guess
	  x0.kpl        = unifrnd(.06,.24);
	  x0.kve_ve      = unifrnd(.02,10);
	  x0.T1P        = unifrnd(20 ,40);
	  x0.T1L        = unifrnd(15 ,35);
	  x0.t0         = unifrnd(0  ,8);
	  % remaining params
	  VIF_scale = VIF_scale_ref;
	case(6) 
	  kpl = optimvar('kpl','LowerBound',0);
	  kve_ve = optimvar('kve_ve','LowerBound',0); 
	  T1P = optimvar('T1P','LowerBound',0); 
	  T1L = optimvar('T1L','LowerBound',0); 
	  t0 = optimvar('t0','LowerBound',0); 
	  VIF_scale = optimvar('VIF_scale','LowerBound',0); 
	  % ground truth 
	  xstar.kpl = kpl_ref;
	  xstar.kve_ve = kve_ve_ref;
	  xstar.T1P = T1P_ref;
	  xstar.T1L = T1L_ref;
	  xstar.t0 = t0_ref;
	  xstar.VIF_scale = VIF_scale_ref;
	  % random initial guess
	  x0.kpl        = unifrnd(.06,.24);
	  x0.kve_ve      = unifrnd(.02,10);
	  x0.T1P        = unifrnd(20 ,40);
	  x0.T1L        = unifrnd(15 ,35);
	  x0.t0         = unifrnd(0  ,8);
	  x0.VIF_scale  = unifrnd(10, 1000);
end

% create state variable
state_var = optimexpr([2, Nscans]);
aux_var = optimexpr([2, Nscans]);

% compute Mxy at first scan
i = 1;
aux_var(:, i) = M0;
VIF_i = [fcn2optimexpr(@(bolusstart)VIF_scale*gampdf(params.TRList(i) - bolusstart, gammaPdfA, gammaPdfB), t0); 0];
state_var(:, i) = sin(params.FaList(:,i)) .* (ve_ref * aux_var(:,i) + (1 - ve_ref) * VIF_i);

for i = 2:Nscans
	TR_i = params.TRList(i) - params.TRList(i-1);
	nsubstep = 5;
	deltat = TR_i / nsubstep;
	
	integratedt = params.TRList(i-1) + deltat * [0.5:1:nsubstep];
	%integratedt = [params.TRList(i-1):deltat:params.TRList(i)] + deltat/2  ;

	% manual define gamma fn
	%t0_diff = (integratedt(1:nsubstep)'-t0) ./ gammaPdfB;
  %integrand = VIF_scale * power(t0_diff, gammaPdfA - 1) .* exp(-t0_diff) ...
  %							./ (gammaPdfB * gamma(gammaPdfA));
  % matlab gamma fn
 	integrand = fcn2optimexpr(@(bolusstart2)VIF_scale*gampdf(integratedt(1:nsubstep)' - bolusstart2, gammaPdfA,gammaPdfB), t0);
  
  % exp(A) where A = TR_i * [-1/T1P - kpl - kve_ve, 0; kpl, -1/T1L]
  A11 = -kpl - kve_ve - 1/T1P;
  A22 = - 1/T1L;% - kve_ve;
  expATR = [ exp(TR_i*A11), 0; ...
  					 kpl*(exp(TR_i*A11) - exp(TR_i*A22)) / (A11 - A22), ...
  					 exp(TR_i*A22)];
  
  %t_seq = TR_i - deltat*[.5:1:nsubstep];
  t_seq = params.TRList(i) - integratedt;
  aifterm_1 = (kve_ve * deltat) * (exp(A11*t_seq) * integrand); 
  aifterm_2 = (kve_ve * deltat) * (kpl * (exp(A11*t_seq) - exp(A22*t_seq))/(A11 - A22)) * integrand;
  aifterm = [aifterm_1; aifterm_2];

	aux_var(:,i) = expATR * aux_var(:,i-1) + aifterm;
	% compute Mxy
	VIF_i = [fcn2optimexpr(@(bolusstart3)VIF_scale*gampdf(params.TRList(i) - bolusstart3, gammaPdfA, gammaPdfB), t0); 0];
	state_var(:,i) = sin(params.FaList(:,i)) .* (ve_ref * aux_var(:,i) + (1 - ve_ref) * VIF_i);
	% compute Mz for next time step
	aux_var(:,i) = cos(params.FaList(:,i)).* aux_var(:,i);
end
truth_state = evaluate(state_var, xstar);
truth_aux = evaluate(aux_var, xstar);

% plotting
figure(1)
plot_voxels = [Voxel_ID]
N_plot_voxels = length(plot_voxels)
for i = 1:N_plot_voxels
	pix_i = plot_voxels(i)
	subplot(N_plot_voxels, 2, 2*(i-1) + 1)
	%plot(truth_state(1, :), 'r', 'DisplayName', 'true') 
	%hold on
	pix_i_str = sprintf('voxel %d', pix_i)
	plot(Mxy_data_voxels(1, :, pix_i, 1), 'b', 'DisplayName', pix_i_str)
	hold off
	
	subplot(N_plot_voxels, 2, 2*(i-1) + 2)
	%plot(truth_state(2, :), 'r', 'DisplayName', 'true') 
	%hold on
	plot(Mxy_data_voxels(2, :, pix_i, 1), 'b', 'DisplayName', pix_i_str)
	hold off
end
drawnow
% subplot(2,1,1)
% plot(truth_state(1, :), 'r', 'DisplayName', 'true') 
% hold on
% plot(Mxy_truth_voxels(1, :, 1), 'b', 'DisplayName', 'voxel 1') 
% hold on
% plot(Mxy_truth_voxels(1, :, 11), 'b--', 'DisplayName', 'voxel 11') 
% hold on
% plot(Mxy_truth_voxels(1, :, 16), 'g', 'DisplayName', 'voxel 16')
% hold on
% plot(Mxy_truth_voxels(1, :, 21), 'g--', 'DisplayName', 'voxel 21') 
% hold on
% plot(Mxy_truth_voxels(1, :, 25), 'o', 'DisplayName', 'voxel 25')
% hold on
% plot(Mxy_truth_voxels(1, :, 30), 'o--', 'DisplayName', 'voxel 30')  
% hold off

% subplot(2,1,2)
% plot(truth_state(2, :), 'r', 'DisplayName', 'true') 
% hold on
% plot(Mxy_truth_voxels(2, :, 1), 'b', 'DisplayName', 'voxel 1') 
% hold on
% plot(Mxy_truth_voxels(2, :, 11), 'b--', 'DisplayName', 'voxel 11') 
% hold on
% plot(Mxy_truth_voxels(2, :, 16), 'g', 'DisplayName', 'voxel 16')
% hold on
% plot(Mxy_truth_voxels(2, :, 21), 'g--', 'DisplayName', 'voxel 21') 
% hold on
% plot(Mxy_truth_voxels(2, :, 25), 'o', 'DisplayName', 'voxel 25')
% hold on
% plot(Mxy_truth_voxels(2, :, 30), 'o--', 'DisplayName', 'voxel 30')
% hold off
% drawnow

% loop over voxels
start_voxel = Voxel_ID;
end_voxel = Voxel_ID;
for id_voxel = start_voxel:end_voxel 
	% loop over samples of noisy data and solve inverse problem
	start_trial = 1;
	end_trial = NumTrials+1;
	if RecTrue == 1
		end_trial = 1;
	%else
	%	start_trial = 2;
	end
	for i_trial = start_trial:end_trial
		fprintf('voxel = %d, Sample = %d \n', id_voxel, i_trial);
		
		% get noise amount
		Mxy_std_i = Mxy_std_vec(:, id_voxel);
		
		% get noisy data
		disp('Generating noisy data')
		% i_trial = 1 is noiseless
		if i_trial > 1
			for j_spc=1:2
				Mxy_data_voxels(j_spc, :, id_voxel, i_trial) = Mxy_data_voxels(j_spc, :, id_voxel, 1) + Mxy_std_i(j_spc) * randn( size( Mxy_data_voxels(j_spc, :, id_voxel, 1) ) );
			end
		end
		
		% cost function
		disp('Defining cost function')
		cost_fn = sum((state_var(1, :) - Mxy_data_voxels(1, :, id_voxel, i_trial)).^2) + sum((state_var(2, :) - Mxy_data_voxels(2, :, id_voxel, i_trial)).^2);

		% create optimization problem
		disp('Create optimization problem')
		opt_prob = optimproblem('Objective', cost_fn);
		% show(opt_prob); % to view
		% problem = prob2struct(opt_prob, 'ObjectiveFunctionName', 'generatedObjective');

		% options
		opt_params = optimoptions(@lsqnonlin, 'Display', 'iter-detailed')%, ...
																% 'FunctionTolerance', 1.e-9, ...
																% 'MaxFunctionEvaluations', 5000, ...
																% 'MaxIterations', 500);

		% initial state
		init_state = evaluate(state_var, x0);
		init_aux = evaluate(aux_var, x0);

		% solve
		disp('Solving ...')
		[popt, fval, exitflag, output] = solve(opt_prob, x0, ...
																'Options', opt_params, ...
																'solver', 'lsqnonlin', ...
																'ObjectiveDerivative', 'finite-differences');

		% store results 
		switch(NumberRecParams)
		case(1)
			rec_params(i_trial, 1, id_voxel) = popt.kpl;
		case(2)
			rec_params(i_trial, 1, id_voxel) = popt.kpl;
			rec_params(i_trial, 2, id_voxel) = popt.kve_ve;
		case(3)
			rec_params(i_trial, 1, id_voxel) = popt.kpl;
			rec_params(i_trial, 2, id_voxel) = popt.kve_ve;
			rec_params(i_trial, 3, id_voxel) = popt.t0;
		case(4)
			rec_params(i_trial, 1, id_voxel) = popt.kpl;
			rec_params(i_trial, 2, id_voxel) = popt.kve_ve;
			rec_params(i_trial, 3, id_voxel) = popt.T1P;
			rec_params(i_trial, 4, id_voxel) = popt.T1L;
		case(5)
			rec_params(i_trial, 1, id_voxel) = popt.kpl;
			rec_params(i_trial, 2, id_voxel) = popt.kve_ve;
			rec_params(i_trial, 3, id_voxel) = popt.T1P;
			rec_params(i_trial, 4, id_voxel) = popt.T1L;
			rec_params(i_trial, 5, id_voxel) = popt.t0;
		case(6)
			rec_params(i_trial, 1, id_voxel) = popt.kpl;
			rec_params(i_trial, 2, id_voxel) = popt.kve_ve;
			rec_params(i_trial, 3, id_voxel) = popt.T1P;
			rec_params(i_trial, 4, id_voxel) = popt.T1L;
			rec_params(i_trial, 5, id_voxel) = popt.t0;
			rec_params(i_trial, 6, id_voxel) = popt.VIF_scale;
		end

		disp('rec_params')
		disp(rec_params(i_trial, :, id_voxel))
	end
end 

% display results
params_mean = mean(rec_params(2:end, :, :));
params_std = std(rec_params(2:end, :, :));
params_var = var(rec_params(2:end, :, :));

params_true = rec_params(1, :, :);

params_mean_kpl = zeros(Nroivoxel,1);
params_std_kpl = zeros(Nroivoxel,1);
for i = 1:Nroivoxel
	params_mean_kpl(i) = mean(rec_params(2:end, 1, i));
	params_std_kpl(i) = std(rec_params(2:end, 1, i));
end

disp('Mean and std of recovered parameters')
for id_voxel = start_voxel:end_voxel
	fprintf('voxel = %d\n', id_voxel);
	disp(params_mean_kpl(id_voxel))
	disp(params_std_kpl(id_voxel))
end

disp('True recovered parameters')
for id_voxel = start_voxel:end_voxel
	fprintf('voxel = %d\n', id_voxel);
	disp(params_true(1, :, id_voxel))
end


% save data
w = whos;
for a = 1:length(w) 
  str.(w(a).name) = eval(w(a).name); 
end
save(fname_data, 'str');  
save(fname_data, '-struct', 'str'); 

