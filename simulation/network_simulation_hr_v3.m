%% Simulate an integrate-and-fire neuron
clear;clc
% Parameters...
sr = 10;      % sampling rate (kHz)
tau     = 20;    % membrane time constant (ms)
tau_AHP = 100;      % after hyperpolarization time constant(ms)
R     = 40;    % membrane resistance (MOhms)
E_L   = -65;   % resting potential for the leak (mV)
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
n = 300; % number of neuron
mu_sp = exp(-.9 + .2*randn(1,n)); % .8 .2
sd = .4;
%%
%common input
I_common = pinknoise_filtered(timestep,Fs,1);
srange = 50*sr;
shift = round(rand(1,n)*srange);

w_com = 0.5;

% random location
for i = 1:n
    xx(i) = 0+1000*rand(1);
    yy(i) = 0+1000*rand(1);
end
for i = 1:n
    for j = 1:n
        dis(i,j) = sqrt((xx(i)-xx(j))^2 + (yy(i)-yy(j))^2);
    end
end
dis_v = dis;
dis_v(dis_v==0) = [];
dis50 = prctile(dis_v,50);

% random v
velocity = rand(n,1)*150 + 60;
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

% figure,
% histogram(rand_logn_exc,'Normalization','probability')
% hold on
% histogram(-rand_logn_inh,'Normalization','probability')
%% random weight assignment
w_syn = zeros(n,n);
numsyn = 0;
tau_syn = 1*sr; % ms
t_syn = 10*sr; % ms
I_syn = zeros(n,n,t_syn);
sign = zeros(1,n);
ii = randperm(n,n*.8); % 80%exc, 20%inh
sign(ii) = 1;
for i = 1:n
    
    if sign(i) == 1
        syni = randperm(n,n*.1); %connection prob = 0.1*50% = 5%
    elseif sign(i) == 0
        sign(i) = -1;
        syni = randperm(n,n*.4); %connection prob = 0.4*50% = 5%
    end
    
    for j = 1:n
       
        if dis(i,j) > dis50 || i == j || sum(ismember(j,syni))==0
            continue
        end
                
        numsyn = numsyn +1;
        if sign(i) > 0
            w_syn(i,j) = rand_logn_exc(numsyn);
        else
            w_syn(i,j) = -rand_logn_inh(numsyn);
        end
        
        I_syn(i,j,:) = w_syn(i,j) * (tau_syn:tau_syn+t_syn-1)/tau_syn.*exp(1-(tau_syn:tau_syn+t_syn-1)/tau_syn); %alpha function, dt = 0
        
    end
end

syn_den = sum(sum(w_syn ~= 0))/(n*(n-1))

fprintf('Basic info finished')

%% simulation
clear p;
dt = 1/sr;
V = zeros(n,1);
Ca = zeros(n,1);
V_temp = -65.001*ones(n,1);
Ca_temp = zeros(n,1);
dt_true = zeros(n,n);
%dt_max = ceil(dis50/min(velocity)*sr);
dt_max = ceil(max(dis(:))/min(velocity)*sr);
spk = {};
k = zeros(n,1);
I_cumsyn = zeros(n,t_syn+dt_max);
I_input = zeros(n,1);
I_sp = zeros(n,Fs*10);

for i = 1:n
    post{i} = find(w_syn(i,:)~=0);
end
for t = 1:timestep
    if mod(t,Fs) == 0
        fprintf('Time %d s...',t/(Fs))
        fprintf('\n')
    end
    % spontaneous firing I + common input (generate this every 10 secs)
    if mod(t,Fs*10) ==1
        for i = 1:n
            I_common_s = circshift(I_common,shift(i)) ;
            I_sp(i,:) = mu_sp(i) + sd*( (1-w_com)*pinknoise_filtered(Fs*10,Fs,1) + w_com*I_common_s(t:t+Fs*10-1) );
        end
    end
    
    if t == 1
        continue;
    end
    % I_input at time t for all the neurons
    % update cumulative synaptic input
    tt = mod(t,Fs*10) + Fs*10*(mod(t,Fs*10)==0);
    for i = 1:n
        I_input(i) = I_sp(i,tt) + I_cumsyn(i,1);
        I_cumsyn(i,:) = [I_cumsyn(i,2:end),0];
    end
    
    for i = 1:n
        
        V(i) = V_temp(i) + ...
            dt/tau*(-(V_temp(i)-E_L)- Ca_temp(i)*g_AHP*(V_temp(i)-V_AHP)*R+I_input(i)*R);
        
        
        % after hyperpolarization current
        Ca(i) = Ca_temp(i) - Ca_temp(i)/tau_AHP*dt;
        
        % Fire...
        if V(i)>Vthres
            k(i) = k(i)+1;
            spk{i}(k(i)) = t;
            % reset
            V(i) = Vreset;
            Ca(i) = Ca_temp(i)+alpha;
            % send synaptic signal
            for j = post{i}
                
                latency = round(dis(i,j)/velocity(i)*sr);
                latency = max(latency,1);
                dt_true(i,j) = latency;
                if t+latency+t_syn-1 > timestep
                    continue
                end
                I_cumsyn(j,(latency):(latency+t_syn-1)) = I_cumsyn(j,(latency):(latency+t_syn-1)) + squeeze(I_syn(i,j,:))';
                
            end
        end
        Ca_temp(i) = Ca(i);
        V_temp(i) = V(i);
        
    end
end
fprintf('Neurons finished!')
%%
data.mu_sp = mu_sp;
data.distance = dis;
data.xx = xx;
data.yy = yy;
data.velocity = velocity;
data.frate = k/(timestep/Fs);
for i = 1:n
    data.spk{i} = spk{i}/sr; %unit: ms
end
data.syn.w_syn = w_syn;
data.syn.syn_den = syn_den;
data.syn.dt_true = dt_true/sr;

save('simulated_data_sim1.mat','data')
fprintf('Neurons saved!')
%% plot CCG
pre = 1;post = 2;
[ccg,deltaT] = corr_fast_v3(data.spk{pre},data.spk{post},-25,25,101); 
ccg = ccg(1:101);
figure,
bar(linspace(-25,25,101),ccg)
xlabel('Interval [ms]')
