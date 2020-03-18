

x = -100:100;
y = x*0;
tau=10;
y = x*0;
y(x>0) = exp(-x(x>0)/tau)-exp(-x(x>0)/(tau/10));
y=y/max(y);
y = y+randn(size(y))/10;
plot(x,y)

yhat = fitXcorrSpline(x',y',10);
hold on
plot(x,yhat,'r')
hold off