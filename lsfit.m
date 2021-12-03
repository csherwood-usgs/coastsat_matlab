function [a,b,r2,sa,sb,hdot] = lsfit(x,y,iplot,t)
% LSFIT  Least-squares fit
%
% [a,b, r2, sa, sb, hdot ] = lsfit(x,y,[iplot])
%
% Least-squares fit of
%   y = a + b * x
% x and y must be column vectors of same length.
% I think sa and sb are standard error of offset and slope.
%
% Makes a plot if iplot == 1, No plot if iplot is absent
% If provided, sd is plotted.
% r2 calculated according to Sachs(1982) p. 389
% hdot is the handle for the data dots in the plot

%    Chris Sherwood, UW and USGS
%    Last revised March 19, 1999

echo off
format compact
if(nargin<3),iplot=0;,end;

% matlab ls method
% x = [x ones(x) ];
% b = x\y

% or, equivalently
% b = (x'*x)\(x'*y)
% b = polyfit( log( z ./ za ), y, 1 )
% b = inv(x'*x)*(x'*y)
% predictions
% yhat = x*b;
% or
% yhat2 = polyval( b2, log( z ./ za ) );

% or, longhand
n = max(size(x));
sumx = sum(x);
sumy = sum(y);
sumx2 = sum(x .* x);
sumy2 = sum(y .* y);
sumxy = sum(x .* y);

% (see Taylor, p. 157)
% del = n*sumx2 - sumx*sumx;
% a = (1 ./del)*(sumx2*sumy - sumx*sumxy)
% b = (1 ./del)*(n*sumxy-sumx*sumy)
% yhat = [x ones(x)]*[ b; a ];
% varY = ( 1 ./(n-2) )*sum( (y-yhat).^2 );
% uncertainty in a = b(2)
% varA = varY * sumx2/del;
% uncertainty in b = b(1)
% varB = n*varY/del;
% sdp = sqrt( varB );
% r2  = sum( (yhat-mean(y)).^2 ) / sum( (y-mean(y)).^2 )

% (Sachs, p. 417)
Qx = sumx2-sumx^2/n;
Qy = sumy2-sumy^2/n;
Qxy = sumxy-(1/n)*sumx*sumy;

xbar = sumx/n;
ybar = sumy/n;

b = Qxy/Qx;
% or b = (sumxy-(sumx*sumy)/n)/(sumx2-sumx^2/n)
a = (sumy-b*sumx)/n;

sx = sqrt(Qx/(n-1));
sy = sqrt(Qy/(n-1));
r = Qxy/sqrt(Qx*Qy);
r2 = r*r;
Qydotx = Qy-b*Qxy;
sydotx = sqrt( Qydotx/(n-2) );
% or sydotx = sy*sqrt( (1-r2)*(n-1)/(n-2) );

% s.d. of intercept and slope
sb = sydotx/sqrt(Qx);
sa = sydotx*(sqrt(1/n+xbar^2/Qx));
%fprintf(1,'Slope: %f +- %f\nIntercept: %f +- %f\n',b,sb,a,sa);

t = 2.132; %n-2=4, two sided, alpha = .9

xline = linspace(min(x),max(x),20);
yline = a+b*xline;
syhatdot = sydotx*sqrt(1+1/n+((xline-xbar).^2)./Qx);
yplus = yline+syhatdot;
yminus = yline-syhatdot;
yhat = a+b*x;
residuals = (y - yhat)';

hdot = nan;
if(iplot==1),
hold on
h=plot(xline,yline,'--k','linewidth',3);
set(h,'color',[.4 .4 .4])
h=plot(xline,yplus,':k','linewidth',2);
set(h,'color',[.4 .4 .4])
h=plot(xline,yminus,':k','linewidth',2);
set(h,'color',[.4 .4 .4])
hdot=plot( x,y,'ok','markersize',4,'markerfacecolor',[.8 .1 .1]);

end






