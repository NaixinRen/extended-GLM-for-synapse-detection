function plotCCG(pre,post,CCG,t,model_fits,results)

NeuronType = results.neurontype;
if numel(fieldnames(results)) ==3
    CnxLabel = results.cnx_label;
end
bar(t,CCG{pre,post},'k')
hold on
p1 = plot(t,model_fits(pre).glm.yhat(post,:),'Color',[0.9290    0.6940    0.1250],'LineWidth',2);
hold on
if NeuronType(pre) > 0
    p2 = plot(t,model_fits(pre).exc.yhat(post,:),'Color',[0.8500    0.3250    0.0980],'LineWidth',2);
    legend([p1 p2],{'Slow Model','Exc Full Model'})
elseif NeuronType(pre) < 0
    p2 = plot(t,model_fits(pre).inh.yhat(post,:),'Color',[     0    0.4470    0.7410],'LineWidth',2);
    legend([p1 p2],{'Slow Model','Inh Full Model'})
end

if numel(fieldnames(results)) ==3
    title([num2str(pre) ' --> ' num2str(post) ', Putaitve Connection: ' num2str(CnxLabel(pre,post))])
else
    title([num2str(pre) ' --> ' num2str(post)])
end
hold off
