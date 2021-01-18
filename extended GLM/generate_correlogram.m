function [CCG, t, distance, ignore] = generate_correlogram(spikes,sr,location,hyperparameter,ignore_index)


% This function generates cross-correlograms from spike trains

% input arguments:
% -spikes: spike times, in MILLISECONDS
% -sr: sampling rate of the spike data (per ms) 
% -location: 2D location of the neurons(soma)
% -interval: interval of the correlogram. eg. interval = 50 means the interval of the correlogram is [-25,25] ms
% -ignore_index: 1 or 0, when ignore_index == 1,ignore low firing rate neurons (# of spikes <1000), sparse correlograms 
% (# of spikes in a correlogram <50), and neuron pairs with possbile spike sorting problems. Defalt value is 1.


if nargin < 4
    error('Not enough arguments')
elseif nargin == 4
    ignore_index = 1;
end

NN = length(spikes); 

distance = zeros(NN,NN);
for i = 1:NN
    for j = 1:NN
        distance(i,j) = sqrt((location.x(i)-location.x(j))^2 + (location.y(i)-location.y(j))^2);
    end
end

Tlist=cell(0);
for i=1:NN
    Tlist{i}=spikes{i};
    %add sub-resolution noise
    Tlist{i}=Tlist{i}+(rand(size(Tlist{i}))-.5)*(1/sr);
end

interval = hyperparameter.interval;
binsize = hyperparameter.binsize;
t = linspace(-interval/2,interval/2,interval/binsize+2);
t = t+mean(diff(t))/2;
t = t(1:interval/binsize+1);

CCG = cell(0);
spikes_pre = zeros(NN,NN);
for i=1:NN
    if length(Tlist{i}) <= 1000 && ignore_index ==1
        continue
    end
    for j=(i+1):NN
        if length(Tlist{j}) <= 1000 && ignore_index ==1
            continue
        end
        [tmp,~] = corr_fast_v3(Tlist{i},Tlist{j},-interval/2,interval/2,interval/binsize+2);
        CCG{i,j} = tmp(1:interval/binsize+1)'; % last bin of histc contains data == the last bin
        CCG{j,i} = flipud(tmp(1:interval/binsize+1))';
    end
end

ignore = eye(NN,NN);
if ignore_index ==1
ignore(cellfun(@sum,CCG) < 50) = 1;
% find neuron pairs with possible spike sorting problems
[ignore_spkst,SortingProbIndex] = spkstprob(CCG,interval,binsize);
ignore(ignore_spkst==1) = 1;
end
