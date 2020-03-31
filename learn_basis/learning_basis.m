
function [X] = learning_basis(CCG,ignore)

addpath(genpath('fastBSpline\'));
addpath(genpath('minFunc\'));
%%

n = size(CCG,2);
ijlut = zeros(n*(n-1)/2,2);
k =0;
for i =1:n-1
    for j = i+1:n
       
        if isempty(CCG{i,j}) || ignore(i,j) == 1
            continue
        end
         k = k+1;
        y(k,:) = CCG{i,j};
        ijlut(k,1) = i;
        ijlut(k,2) = j;
    end
end

%%

[n,m]=size(y);
k=6; % rank to approximate

% set up bases
pb = 16;
bbas = getCubicBSplineBasis(linspace(0,1,m),pb,0);
bbas = bbas(:,2:end);

% random initialization
B = randn(k,pb);

%%

X = B*bbas';
A = zeros(n,k+1);
% w_center = ones(101,1);
% w_center(49:52,1) = [0;0;0;0];
% w_center = w_center';
for i=1:n
    idx = isfinite(y(i,:));
    [A(i,:)] = glmfit(X(:,idx)',y(i,idx)','poisson'); % ,'weights',w_center);
end
mu=A(:,1); A=A(:,2:end);
%%
options=[];
options.method = 'cg';
options.Display = 'off';
[x,f,exitflag,output] = minFunc(@lossLowRank,[mu; A(:); B(:)],options,y,k,bbas');
mu = x(1:n);
A = reshape(x((n+1):(n*k+n)),n,k);
B = reshape(x((n*k+n+1):end),k,pb);
X = B*bbas';
%X = X';


