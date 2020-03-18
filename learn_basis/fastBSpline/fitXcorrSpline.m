

function yhat = fitXcorrSpline(x,y,nknots)

% Smooth using B-Splines for acorr...
x0 = [-1 linspace(0,1,length(x)) 2];
knots = sort([-1 -logspace(-3,0,nknots)/2 logspace(-3,0,nknots)/2 1])+0.5;

s = fastBSpline.pspline(knots,3,x0,[0; y; 0],2);
yhat = s.evalAt(x0(2:(end-1)));