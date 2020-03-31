function b1 = random_parameter_s1(y,XX,t,b_glm,distance,eta,tau0,stage)

% t = linspace(-25,25,102);
% t = t+mean(diff(t))/2;
% t = t(1:101);
v0 = [0;0];


options=[];
options.method = 'cg';
options.MaxIter = 50;
options.Display = 'off';
f=Inf;

switch stage
    
    case 1
        
        for rr=1:50
            rlat = log(gamrnd(2,2)/2+.01);
            if mod(rr,2)==0
                b0 = [log(nanmean(y)), b_glm(2:end)*0,log(abs(randn(1))), rlat, log(tau0)];
            else
                b0 = [b_glm(1), b_glm(2:end)+randn(1,size(XX,1)-1)/5, log(abs(randn(1))), rlat, log(tau0)];
            end
            [brr,frr] = minFunc(@loss_excalpha,b0',options,XX',y',t',v0,distance,eta,tau0);
            if frr<f
                b1=brr;
                f=frr;
            end
        end
        
    case -1
        
         for rr=1:50
            rlat = log(gamrnd(2,2)/2+.01);
            if mod(rr,2)==0
                b0 = [log(nanmean(y)), b_glm(2:end)*0,log(abs(randn(1))), rlat, log(tau0)];
            else
                b0 = [b_glm(1), b_glm(2:end)+randn(1,size(XX,1)-1)/5, log(abs(randn(1))), rlat, log(tau0)];
            end
            [brr,frr] = minFunc(@loss_inhalpha,b0',options,XX',y',t',v0,distance,eta,tau0);
            if frr<f
                b1=brr;
                f=frr;
            end
        end
end