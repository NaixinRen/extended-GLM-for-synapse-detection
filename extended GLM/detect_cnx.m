function results = detect_cnx(model_results,ignore,threshold)

% This function returns the final results ("results") based on the model fits. 
%
% results is a structure with 3 fields:
%
% -cnx_label: the detected putative connections based on the model results
%   and threshold. 1 = putative exc cnx, -1 = putative inh cnx.
%
% -neurontype: the presynaptic neuron type based on the model results. 1 =
%   putative exc neuron, -1 = putative inh neuron.
%
% -llr_matrix: the log likelihood ratio (llr) between full model and slow model
%   for each neuron pair. This is the score that used to detect synaptic cnx.
%   If llr > threshold, the neuron pair is a putative connection.
%
% If the user doesn't specify a threshold, the function will only return
% neurontype and llr_matrix


if nargin < 2
    error('Not enough arguments')
elseif nargin == 2
    threshold = Inf;
end


NN = size(model_results,2);
llr_max = nan(NN,1);
neurontype = nan(NN,1);
llr_matrix = nan(NN,size(ignore,2));
cnx_label = nan(NN,size(ignore,2));

for pre = 1:NN
    
    if sum(1-ignore(pre,:)) == 0
        continue
    end
    
    exc_llr = model_results(pre).exc.llh-model_results(pre).glm.llh;
    inh_llr = model_results(pre).inh.llh-model_results(pre).glm.llh;
    llr = exc_llr-inh_llr;
    
    [~,ii] = nanmax(abs(llr));
    llr_max(pre,1) = llr(ii);
    
    if llr_max(pre,1) > 0
        neurontype(pre) = 1;
        llr_matrix(pre,:) = exc_llr;
        cnx_label(pre,:) = llr_matrix(pre,:) > threshold;
    elseif llr_max(pre,1) < 0
        neurontype(pre) = -1;
        llr_matrix(pre,:) = inh_llr;
        cnx_label(pre,:) = ( llr_matrix(pre,:) > threshold )*(-1);
    end
    
end
if nargin == 3
    results.cnx_label = cnx_label;
end
results.neurontype = neurontype;
results.llr_matrix = llr_matrix;


