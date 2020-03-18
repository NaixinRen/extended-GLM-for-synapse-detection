function plotCCG(pre,post,CCG,model_fits,results)
NeuronType = results.neurontype;
CnxLabel = results.cnx_label;
t = linspace(-25,25,101);
bar(t,CCG{pre,post},'k')
hold on
p1 = plot(t,model_fits(pre).glm.yhat(post,:),'Color',[0.9290    0.6940    0.1250],'LineWidth',2);
hold on
if NeuronType(pre) > 0
   p2 = plot(t,model_fits(pre).exc.yhat(post,:),'Color',[0.8500    0.3250    0.0980],'LineWidth',2);
elseif NeuronType(pre) < 0
    p2 = plot(t,model_fits(pre).inh.yhat(post,:),'Color',[     0    0.4470    0.7410],'LineWidth',2);
end
legend([p1 p2],{'Slow Model','Full Model'})
title([num2str(pre) ' --> ' num2str(post) ', Connection label: ' num2str(CnxLabel(pre,post))])
hold off
