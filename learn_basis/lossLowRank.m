
function [f,dx] = lossLowRank(b,y,k,bas,offset)

if nargin<5, offset=0; end

[n,m]=size(y);
pb=size(bas,1);

mu = b(1:n);
A = reshape(b((n+1):(n*k+n)),n,k);
B = reshape(b((n*k+n+1):end),k,pb);

xb = mu + A*B*bas + offset;
lam = exp(xb);
w_center = ones(101,1);
w_center(49:52,1) = [0;0;0;0];
w_center = w_center';
f = -nansum(nansum((y.*xb - lam).*w_center));
lam_err = (lam-y).*w_center;
lam_err(~isfinite(lam_err))=0;
dmu = sum(lam_err,2);
dA = (lam_err)*(B*bas)';
dB = A'*(lam_err)*bas';
dx = [dmu(:); dA(:); dB(:)];

