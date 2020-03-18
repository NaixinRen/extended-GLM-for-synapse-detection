function results = detect_cnx(model_results,threshold,ignore)

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

results.cnx_label = cnx_label;
results.neurontype = neurontype;
results.llr_matrix = llr_matrix;


