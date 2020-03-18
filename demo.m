clear; clc;
addpath(genpath('extended GLM\'));
addpath(genpath('learn_basis\'));
load('simulated data_50 neurons.mat')
%% parameter setting
spikes = data.spk;
sr = 10; % sampling rate per ms (kHz)
location.x = data.xx;
location.y = data.yy;
ignore_index = 1; % ignore low firing neurons and neuron pairs with potential spike sorting problem

[CCG, distance, ignore] = generate_correlogram(spikes,sr,location,ignore_index); % convert spike data into cross-correlograms

X = learning_basis(CCG,ignore); % learn basis function
%% model fitting and results
NN = size(CCG,1);
% for parallel computing
% poolobj = gcp('nocreate');
% delete(poolobj);
% parpool('local',2);
for pre = 1:NN % use parfor for parallel computing
    fprintf('neuron %i ',pre)
    model_fits(pre) = extendedGLM(CCG(pre,:), X, distance(pre,:), ignore(pre,:));   
end

threshold = 5.09;
results = detect_cnx(model_fits,threshold,ignore);
%% plot CCG and fitting results

pre = 50;post = 11;
figure,
plotCCG(pre,post,CCG,model_fits,results)

