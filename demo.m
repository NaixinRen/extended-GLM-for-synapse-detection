clear; clc;
addpath(genpath('extended GLM\'));
addpath(genpath('learn_basis\'));
load('simulated data_50 neurons.mat') % load simulated data
%% parameter setting

spikes = data.spk; %spike trains, in MILLISECONDS
sr = 10; % sampling rate per ms (kHz)
location.x = data.xx; % neuron location (x coordinate)
location.y = data.yy; % neuron location (y coordinate)
ignore_index = 1; % ignore low firing neurons and neuron pairs with potential spike sorting problem
hyperparameter.binsize = .5; % binsize of the CCG (ms)
hyperparameter.interval = 50; % interval of the CCG (ms)
hyperparameter.tau0 = 0.8; %ms
hyperparameter.eta_w = 5;
hyperparameter.eta_tau = 20;
hyperparameter.eta_dt_coeff = 2;

[CCG, t, distance, ignore] = generate_correlogram(spikes,sr,location,hyperparameter,ignore_index); 
X = learning_basis(CCG,ignore); 
%% model fitting and results

NN = size(CCG,1);

% for parallel computing
% poolobj = gcp('nocreate');
% delete(poolobj);
% parpool('local',2);

for pre = 1:NN % use parfor for parallel computing
    fprintf('neuron %i ',pre)
    model_fits(pre) = extendedGLM(CCG(pre,:), X, distance(pre,:),hyperparameter,ignore(pre,:));  
end

threshold = 5.09;
results = detect_cnx(model_fits,ignore,threshold);
%% plot CCG and fitting results

pre = 1;post = 14;
figure(1),
plotCCG(pre,post,CCG,t,model_fits,results)
%% ROC curve for synapse detection(regardless of the sign of connections)

true_label = data.syn.w_syn~=0;
true_label = true_label(eye(NN)~=1);
scores = results.llr_matrix(eye(NN)~=1);
scores(isnan(scores)) = -Inf;
[X,Y,T,AUC] = perfcurve(true_label,scores,1);

%% Compare overall connectivity matrices

figure,
subplot(1,3,1)
plot(X,Y,'LineWidth',2)
ylabel('True positive rate')
xlabel('False positive rate')
title('ROC for synapse detection')
axis square

subplot(1,3,2)
imagesc(sign(data.syn.w_syn))
ylabel('Presynaptic Neuron')
xlabel('Postsynaptic Neuron')
title ('true connections')
axis square
subplot(1,3,3)
imagesc(results.cnx_label)
xlabel('Postsynaptic Neuron')
title ('estimated connections')
axis square
