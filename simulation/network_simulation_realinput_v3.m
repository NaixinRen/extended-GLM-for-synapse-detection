%% Simulate an integrate-and-fire neuron
clear;clc
%for parallel computing
poolobj = gcp('nocreate');
delete(poolobj);
parpool('local',16);
%%
clear;clc;
load('DataSet16.mat')
predata{1} = data;
load('DataSet23.mat')
predata{2} = data;
%% Parameters...
sr = 20;      % sampling rate (kHz)
tau    = 20;    % membrane time constant (ms)
tau_AHP = 100;    % after hyperpolarization time constant(ms)
R      = 40;    % membrane resistance (MOhms)
E_L    = -65;   % resting potential for the leak (mV)
V_AHP = -80;   % hyperpolarization voltage (mV)
alpha = 0.2;    % Calcium influx per spike (miuM)
g_AHP = 0.015;  % Conductivity for AHP (mS/cm^2)
g_AHP = g_AHP*10; % S/m^2
Vthres = -50;   % spike threshold (mV)
Vreset = -65;   % reset potential (mV)
Fs = 1000*sr;      % sampling rate (s)
T = 60*60;   % simulation time (s)
timestep = Fs*T;   % simulation timestep (ms)
I_rheobase = (Vthres - E_L)/R %nA
%% spontaneous firing current (units of nA)...
% number of presynaptic neuron is fixed to be 300
n = 300; % number of postsynaptic neuron
mu_sp = exp(-1 + .2*randn(1,n)); % .8 .2
sd = .4;
%%
% shift in common input
srange = 50*sr;
shift = round(rand(1,n)*srange);

% random location
for i = 1:n
    xx(i) = -1000+2000*rand(1);
    yy(i) = -600+1200*rand(1);
end
for i = 1:300
    prexx(i) = -1000+2000*rand(1);
    preyy(i) = -600+1200*rand(1);
end
for pre = 1:300
    for post = 1:n
        dis(pre,post) = sqrt((prexx(pre)-xx(post))^2 + (preyy(pre)-yy(post))^2);
    end
end
dis_v = dis;
dis_v(dis_v==0) = [];
dis50 = prctile(dis_v(:),50);

% random v
velocity = rand(300,1)*200 + 100;
%% synaptic weight distribution ~ lognormal (generate lognomally distributed random numbers)

% exc
w_syn_min = 0.05; 
w_syn_max = .4;
rand_norm = .5*randn(100000,1)-2.5;%+(log(w_syn_min)+log(w_syn_max))/2;
rand_norm = rand_norm(rand_norm<log(w_syn_max) & rand_norm>log(w_syn_min));
rand_logn_exc = exp(rand_norm);

%inh
w_syn_min = 0.05; 
w_syn_max = .4;
rand_norm = .5*randn(100000,1)-2.5;%+(log(w_syn_min)+log(w_syn_max))/2;
rand_norm = rand_norm(rand_norm<log(w_syn_max) & rand_norm>log(w_syn_min));
rand_logn_inh = exp(rand_norm);

%% random weight assignment
w_syn = zeros(300,n);
numsyn = 0;
tau_syn = 1*sr; % ms
t_syn = 10*sr; % ms
I_syn = zeros(300,n,t_syn);
sign = zeros(1,n);
ii = randperm(300,300*.8); % 80%exc, 20%inh
sign(ii) = 1;
for pre = 1:300
    
    if sign(pre) == 1
        syni = randperm(n,n*.1); %connection prob = 0.1*50%
    elseif sign(pre) == 0
        sign(pre) = -1;
        syni = randperm(n,n*.4); %connection prob = 0.4*50%
    end
    
    for post = 1:n
        if dis(pre,post) > dis50 || sum(ismember(post,syni))==0
            continue
        end
        
        
        numsyn = numsyn +1;
        if sign(pre) > 0
            w_syn(pre,post) = rand_logn_exc(numsyn);
        else
            w_syn(pre,post) = -rand_logn_inh(numsyn);
        end
        %I_syn(pre,post,:) = w_syn(pre,post) * (tau_syn:tau_syn+t_syn-1)/tau_syn.*exp(1-(tau_syn:tau_syn+t_syn-1)/tau_syn); %alpha function, dt = 0
        
    end
end


%presnaptic neuron spike train
prespk = [predata{1}.spikes;predata{2}.spikes];
frate = cellfun(@length,prespk);
[~,index] = sort(frate,'descend');
spk=cell(0);
for pre = 1:300
    i = index(pre);
    spk{pre} = floor(prespk{i}*sr);
end
fprintf('Basic info finished')
%% simulation
clear p;
dt = 1/sr;


dt_true = zeros(300,n);
dt_max = ceil(max(dis(:))/min(velocity)*sr);
postspk = {};

parfor post = 1:n
    V_temp = -65.001;
    Ca_temp = 0;
    dtpost = nan(300,1);
    fprintf("neuron %i \n",post)
    
    
    prelist = find(w_syn(:,post)~=0);
    % common input
    [~,II] = sort(dis(:,post));
    spk_bino = zeros(10,timestep);
    for i = 1:10
        spk_bino(i, spk{II(i)}) = 1;
    end
    spk_common = mean(spk_bino);
    spk_common =  filter(normpdf(linspace(-1,1,25*sr)),1,spk_common);
    I_input = mu_sp(post) + sd*( pinknoise_filtered(timestep,Fs,1) + circshift(5*spk_common,shift(post)) );
    fprintf("initialization...\n")
    if isempty(prelist) == 0
        for i = 1:length(prelist)
            
            pre = prelist(i);
            fprintf("preneuron %i \n",pre)
            syndt = round(dis(pre,post)/velocity(pre)*sr);
            syndt = max(syndt,1);
            dtpost(pre) = syndt/sr;
            I_preinput = zeros(1,timestep);
            for j = 1:length(spk{pre})
                if spk{pre}(j)>timestep-t_syn+1
                    continue
                end
                inputtime = spk{pre}(j):spk{pre}(j)+t_syn-1;
                I_preinput(inputtime) = I_preinput(inputtime) + (tau_syn:tau_syn+t_syn-1)/tau_syn.*exp(1-(tau_syn:tau_syn+t_syn-1)/tau_syn);
            end
            I_syn = circshift(I_preinput,syndt);
            I_syn(1:syndt) = 0;
            I_input = I_input + w_syn(pre,post)*I_syn;
        end
    end
    k = 0;
    for t = 2:timestep
                if mod(t,Fs*10) == 0
                    fprintf('Time %d s...',t/(Fs))
                    fprintf('\n')
                end
        
        V = V_temp + ...
            dt/tau*(-(V_temp-E_L)- Ca_temp*g_AHP*(V_temp-V_AHP)*R+I_input(t)*R);
        
        % after hyperpolarization current
        Ca = Ca_temp - Ca_temp/tau_AHP*dt;
        
        % Fire...
        if V>Vthres
            k = k+1;
            postspk{post}(k) = t;
            % reset
            V = Vreset;
            Ca = Ca_temp+alpha;
        end
        Ca_temp = Ca;
        V_temp = V;
        
    end
    dt_true(:,post) = dtpost;
end

fprintf('Neurons finished!')
%%
simdata.mu_sp = mu_sp;
simdata.distance = dis;
simdata.xx = xx;
simdata.yy = yy;
simdata.velocity = velocity;
for i = 1:length(postspk)
    simdata.spk{i} = postspk{i}/sr; % unit: ms
    simdata.frate(i,:) = length(postspk{i})/(timestep/Fs);
end
simdata.prespk = cellfun(@(x) x/sr,spk,'un',0); % unit: ms
simdata.preneuron = index(1:300);
simdata.realdatasource  = 'DataSet16,23 from SSC3';
simdata.sr = sr;
simdata.syn.w_syn = w_syn;
simdata.syn.syn_den = syn_den;
simdata.syn.dt_true = dt_true;

save('simulated_data_realinput.mat','simdata')
fprintf('Neurons saved!')
%% plot CCG
pre = 6;post = 1;
[ccg,deltaT] = corr_fast_v3(simdata.prespk{pre},simdata.spk{post},-25,25,101); 
ccg = ccg(1:101);
figure,
bar(linspace(-25,25,101),ccg,'k')
xlabel('Interval [ms]')