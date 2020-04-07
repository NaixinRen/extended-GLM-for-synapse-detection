function model_fits = extendedGLM(CCG, X, distance,hyperparameter,ignore)

% This function fits all the cross-correlogrms from one presynaptic neuron.
%
% model_fits is a structure with 4 fields:
%
% -glm: the model fitting results of the slow model
% ---yhat: fitted CCG
% ---llh: log likelihood 
% ---b: parameter
%
% -exc: the model fitting results of the excitatory full model (stage 2,
% with latency constraints)
% ---eta: [eta_w eta_dt eta_tau], eta_dt is calculated based on the
% estimation 
% ---v: estimated "conduction velocity", dt = distance/v(1)+v(2)
%
% -inh: the model fitting results of the inhibitory full model (stage 2,
% with latency constraints). The fields have the same meaning with field exc.
%
% -stage1: the model fitting results of the full model (stage 1), with a
% field for exc and inh stage 1 respectively.



%%
interval = hyperparameter.interval;
binsize = hyperparameter.binsize;
tau0 = hyperparameter.tau0;
eta_w = hyperparameter.eta_w;
eta_tau = hyperparameter.eta_tau;
eta_dt_coeff = hyperparameter.eta_dt_coeff;


NN = length(CCG);
t = linspace(-interval/2,interval/2,interval/binsize+2);
t = t+mean(diff(t))/2;
t = t(1:interval/binsize+1);

XX = [X;ones(1,interval/binsize+1)];



eta = [eta_w 0 eta_tau];

%% excitatory model
fprintf('exc model...')
% slow model & full model stage 1

v0 = [0;0];


b_glm = nan(NN,7);
yhat_glm = nan(NN,length(t));
llh_glm = nan(NN,1);

b_s1 = nan(NN,10);
yhat_s1 = nan(NN,length(t));
llh_s1 = nan(NN,1);

for post= 1:NN
    
    if ignore(post)
        continue
    end
    
    y = CCG{post};
    
    % slow model
    b_glm(post,:) = glmfit(XX',y','poisson','constant','off');
    yhat_glm(post,:) = exp(b_glm(post,:)*XX);
    llh_glm(post) = nansum((y.*log( yhat_glm(post,:)+( yhat_glm(post,:)==0))- yhat_glm(post,:)));
    
    % full model Stage 1 (without latency constraints):
   
    b1 = random_parameter_s1(y,XX,t,b_glm(post,:),distance(post),eta,tau0,1); % 50 times random restart  
    
    options=[];
    options.method = 'cg';
    options.MaxIter = 500;
    options.Display = 'off';
    [b,~] = minFunc(@loss_excalpha,b1,options,XX',y',t',v0,distance(post),eta,tau0);
    [f,df,lam] = loss_excalpha(b,XX',y',t',v0,distance(post),eta,tau0);
    b_s1(post,:) = b';
    yhat_s1(post,:) = lam';
    llh_s1(post) = nansum((y.*log(lam'+(lam'==0))-lam'));
    
end

% Estimation of the “conduction velocity”

llr = llh_s1 - llh_glm;

llr(isnan(llr))=-Inf;
[~,ii_l] = sort(llr,'descend');
[~,l_order] = sort(ii_l);
weight = 1./(1+exp(2*(l_order-5)));
llr(llr==-Inf)=nan;


[v1,f1] = minFunc(@loss_velocity,[log(1/200);log(0.1)],options,distance',exp(b_s1(:,size(XX,1)+2)),weight);
[f1,df1,err] = loss_velocity(v1,distance',exp(b_s1(:,size(XX,1)+2)),weight);
[ww,ii]=sort(weight,'descend');
se_velocity = sqrt(nansum(err(ii(1:5)).^2)/5);
v = exp(v1)';
eta_dt = eta_dt_coeff/se_velocity;
eta = [eta_w eta_dt eta_tau];


% full model Stage 2 (with latency constraints)

b_s2 = nan(NN,10);
yhat_s2 = nan(NN,length(t));
llh_s2 = nan(NN,1);

for post= 1:NN
    
    if ignore(post)
        continue
    end
    
    y = CCG{post};
    
    % random restart
    
    b1 = random_parameter_s2(y,XX,t,b_s1(post,:),distance(post),v,eta,tau0,1);    
    
    [b,f] = minFunc(@loss_excalpha,b1,options,XX',y',t',v,distance(post),eta,tau0);
    [f,df,yhat_s2(post,:)] = loss_excalpha(b,XX',y',t',v,distance(post),eta,tau0);
    b_s2(post,:) = b';
    lam = yhat_s2(post,:);
    llh_s2(post,1) = nansum((y.*log(lam+(lam==0))-lam));
end

%% output

model_fits.glm.yhat = yhat_glm;
model_fits.glm.llh = llh_glm;
model_fits.glm.b = b_glm;

model_fits.exc.yhat = yhat_s2;
model_fits.exc.llh = llh_s2;
model_fits.exc.b = b_s2;
model_fits.exc.eta = eta;
model_fits.exc.v = v;

model_fits.inh = [];

model_fits.stage1.exc.yhat = yhat_s1;
model_fits.stage1.exc.llh = llh_s1;
model_fits.stage1.exc.b = b_s1;

%% inhibitory model
fprintf('inh model...\n')
% full model stage 1

b_s1 = nan(NN,10);
yhat_s1 = nan(NN,length(t));
llh_s1 = nan(NN,1);
eta = [eta_w 0 eta_tau];

for post= 1:NN
    
    if ignore(post)
        continue
    end
    
    y = CCG{post};
    
    % full model Stage 1 (without latency constraints):
    
    b1 = random_parameter_s1(y,XX,t,b_glm(post,:),distance(post),eta,tau0,-1); % 50 times random restart  
    
    options=[];
    options.method = 'cg';
    options.MaxIter = 500;
    options.Display = 'off';
    [b,f] = minFunc(@loss_inhalpha,b1,options,XX',y',t',v0,distance(post),eta,tau0);
    [f,df,yhat_s1(post,:)] = loss_inhalpha(b,XX',y',t',v0,distance(post),eta,tau0);
    b_s1(post,:) = b';
    lam = yhat_s1(post,:);
    llh_s1(post) = nansum((y.*log(lam+(lam==0))-lam));
    
end

% Estimation of the “conduction velocity”

llr = llh_s1 - llh_glm;

llr(isnan(llr))=-Inf;
[~,ii_l] = sort(llr,'descend');
[~,l_order] = sort(ii_l);
weight = 1./(1+exp(2*(l_order-5)));
llr(llr==-Inf)=nan;


[v1,f1] = minFunc(@loss_velocity,[log(1/200);log(0.1)],options,distance',exp(b_s1(:,size(XX,1)+2)),weight);
[f1,df1,err] = loss_velocity(v1,distance',exp(b_s1(:,size(XX,1)+2)),weight);
[ww,ii]=sort(weight,'descend');
se_velocity = sqrt(nansum(err(ii(1:5)).^2)/5);
v = exp(v1)';
eta_dt = eta_dt_coeff/se_velocity;
eta = [eta_w eta_dt eta_tau];


% full model Stage 2 (with latency constraints)

b_s2 = nan(NN,10);
yhat_s2 = nan(NN,length(t));
llh_s2 = nan(NN,1);

for post= 1:NN
    
    if ignore(post)
        continue
    end
    
    y = CCG{post};
    
    % random restart
    
    b1 = random_parameter_s2(y,XX,t,b_s1(post,:),distance(post),v,eta,tau0,-1);    
    
    [b,f] = minFunc(@loss_inhalpha,b1,options,XX',y',t',v,distance(post),eta,tau0);
    [f,df,yhat_s2(post,:)] = loss_inhalpha(b,XX',y',t',v,distance(post),eta,tau0);
    b_s2(post,:) = b';
    lam = yhat_s2(post,:);
    llh_s2(post,1) = nansum((y.*log(lam+(lam==0))-lam));
end
%% output

model_fits.inh.yhat = yhat_s2;
model_fits.inh.llh = llh_s2;
model_fits.inh.b = b_s2;
model_fits.inh.eta = eta;
model_fits.inh.v = v;

model_fits.stage1.inh.yhat = yhat_s1;
model_fits.stage1.inh.llh = llh_s1;
model_fits.stage1.inh.b = b_s1;
