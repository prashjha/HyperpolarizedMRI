clear all
close all
clf
clc

addpath ../.

% set seed
rng(2) 

% set filenames for log and optimization results
fdir = './save';
if ~exist(fdir, 'dir')
       mkdir(fdir);
end
fname = sprintf('%s%s', fdir, '/noise_calculation'); % <-- set filename
fname_diary = sprintf('%s%s', fname, '.log');
fname_der = sprintf('%s%s', fname, '.mat');
fname_der_str = sprintf('%s%s', fname, '_struct.mat');
diary(fname_diary) % <-- log terminal output

%% model and variable Setup
model = HPModel();

% get default params (per discussion with Jim, Chris, David)
params = model.defaultParams();
M0 = [0; 0];

%% get reference noise
[t, Mxy, Mz] = model.compile(M0, params);

Mxy_pyr_max = max(Mxy(1,:));
Mxy_lac_max = max(Mxy(2,:));

fprintf('\nMxy_pyr_max = %10.8f, Mxy_lac_max = %10.8f\n', Mxy_pyr_max, Mxy_lac_max)

%% save data
w = whos;
for a = 1:length(w) 
  str.(w(a).name) = eval(w(a).name); 
end
save(fname_der_str, 'str');  
save(fname_der, '-struct', 'str'); 

%% compute noisy data and plot them
snrs = [2, 5, 10, 15, 20];
sigma_refs = Mxy_pyr_max./snrs;
fprintf('sigma_refs = ');
for i=1:length(snrs)
    fprintf('%8.6f', sigma_refs(i));
    if i == length(snrs)
        fprintf('\n\n');
    else
        fprintf(', ');
    end
end 

figure(1)
subplot(2,1,1)
plot(t, Mz(1,:), t, Mz(2,:))
legend('Pyruvate Mz', 'Lactate Mz')

subplot(2,1,2)
plot(t, Mxy(1,:), t, Mxy(2,:))
legend('Pyruvate Mxy', 'Lactate Mxy')
set(findall(gcf,'-property','FontSize'),'FontSize',14)

saveas(gcf,sprintf('%s%s', fdir, '/true.png'));

%% generate noisy data
noisy_G_vec = [];
for snr_i = 1:length(snrs)
    figure(snr_i + 1)   
    suptitle( sprintf( 'Mxy plot of noisy data for SNR = %d', snrs(snr_i) ) )

    sigma_ref = sigma_refs(snr_i);

    % generate 10 samples of noisy data
    for j=1:10
        noise = randn(size(Mxy));
        % noise = noise .* noise;
        data = Mxy(:, :) + sigma_ref * noise;
        subplot(5,2,j)
        plot(t, Mxy(1,:), t, Mxy(2,:))
        hold on
        plot(t, data(1,:), t, data(2,:))
        if j == 1
            leg = legend('Pyr', 'Lac', ...
               'Noisy Pyr', 'Noisy Lac'); %, ...
               %'NumColumns', 2, 'Location', 'east');
            %legend('boxoff')
            leg.ItemTokenSize = [10,10];
        end
        set(findall(gcf,'-property','FontSize'),'FontSize',8)
        saveas(gcf,sprintf('%s%s_snr_%d_noise_%8.6f.png', fdir, ...
                           '/noisy_data', snrs(snr_i), sigma_ref));
    end

end

