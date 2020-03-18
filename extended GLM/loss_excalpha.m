function [f,df,lam] = loss_excalpha(b,X,y,t,v,dis,eta)

p = b((size(X,2)+1):end);
p=exp(p);
b = b(1:size(X,2));

w = p(1);
deltat = p(2); 
tau = p(3);
eta_w = eta(1);
eta_dt = eta(2);
eta_tau = eta(3);

t(t<deltat)=deltat;
syn = w*(t-deltat)/tau.*exp(1-(t-deltat)/tau);
mm = max(abs(syn)) +(syn==0);
syn = syn./mm*abs(w);

lam = exp(X*b + syn);
f = -nansum((y.*log(lam+(lam==0))-lam));
fw = eta_w*(w-0)^2;
fdt = eta_dt*((deltat)-dis*v(1)-v(2))^2; 
ftau = eta_tau*(tau-.8)^2;
f = f+fw+fdt+ftau;

dl = (lam - y);
dl(~isfinite(dl)) = 0;
db = X'*dl;
dw = syn'*dl/w+eta_w*2*(w-0); 
ddeltat = nansum((syn/tau-syn./(t-deltat)).*dl)+eta_dt*2*((deltat)-dis*v(1)-v(2));
dtau = nansum(((-syn/tau+syn.*(t-deltat)/tau.^2)).*dl)+eta_tau*2*(tau-.8);
df =  [db;dw*w;ddeltat*deltat;dtau*tau]; 