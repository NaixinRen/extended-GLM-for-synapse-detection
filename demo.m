clear; clc;
addpath(genpath('extended GLM\'));
addpath(genpath('learn_basis\'));
load('simulated data_20 neurons.mat') % load simulated data
%% parameter setting
spikes = data.spk;
sr = 10; % sampling rate per ms (kHz)
location.x = data.xx;
location.y = data.yy;
hyperparameter.binsize = .5; % binsize of the CCG (ms)
hyperparameter.interval = 50; % interval of the CCG (ms, interval=50 means the interval is [-25,25] ms) 
ignore_index = 1; % ignore low firing neurons and neuron pairs with potential spike sorting problem

[CCG, t, distance, ignore] = generate_correlogram(spikes,sr,location,hyperparameter,ignore_index); % convert spike data into cross-correlograms

X = learning_basis(CCG,ignore); % learn basis function

hyperparameter.tau0 = 0.8; %ms
hyperparameter.eta_w = 5;
hyperparameter.eta_tau = 20;
hyperparameter.eta_dt = 2;

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
%%
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
figure,
plot(X,Y,'LineWidth',2)
ylabel('True positive rate')
xlabel('False positive rate')
title('ROC for synapse detection')
%% Compare overall connectivity matrices

figure,
subplot(1,2,1)
imagesc(data.syn.w_syn)
ylabel('Presynaptic Neuron')
xlabel('Postsynaptic Neuron')
subplot(1,2,2)
imagesc(results.cnx_label)
xlabel('Postsynaptic Neuron')

