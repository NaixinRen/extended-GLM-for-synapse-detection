function b1 = random_parameter_s2(y,XX,b_s1,distance,v,eta,stage)

t = linspace(-25,25,102);
t = t+mean(diff(t))/2;
t = t(1:101);
tau0 = 0.8;

options=[];
options.method = 'cg';
options.MaxIter = 50;
options.Display = 'off';
f=Inf;


dt_st = min(distance*v(1)+v(2),20);
synt = t;
synt(t<dt_st)=dt_st;
syn = stage*(synt-dt_st)/tau0.*exp(1-(synt-dt_st)/tau0);

b_est= glmfit([XX',syn(1:101)'],y','poisson','constant','off');
w_st = abs(b_est(end))+(b_est(end)==0)*.001;

switch stage
    
    case 1
        
        for rr = 1:2
            if mod(rr,2) == 0
                b0 = [b_est(1:end-1);log(w_st); log(dt_st); log(tau0)];
            else
                b0 = b_s1';
            end
            
            [brr,frr] = minFunc(@loss_excalpha,b0,options,XX',y',t',v,distance,eta);
            if frr<f
                b1=brr;
                f=frr;
            end
        end
        
    case -1
        
        for rr = 1:2
            if mod(rr,2) == 0
                b0 = [b_est(1:end-1);log(w_st); log(dt_st); log(tau0)];
            else
                b0 = b_s1';
            end
            
            [brr,frr] = minFunc(@loss_inhalpha,b0,options,XX',y',t',v,distance,eta);
            if frr<f
                b1=brr;
                f=frr;
            end
        end
        
end

