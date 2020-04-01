function h=patchell(A,b,plotType,color)
% [xe,ye]=plotellip(A,xc,linetype,couleur);
%   plots ellipsoid (x-xc)'*A*(x-xc) = 1 in R^2, or
% [xe,ye]=plotellip(A,b,linetype,couleur);
%   plots ellipsoid x'*A*x + b'*x + c = 0 in R^2

if nargin==2
    type1=1;  xc=b;  linetype='-'; color = 'k';
elseif nargin==4
    type1=1;  xc=b;
else
    type1=0;
end

A=.5*(A+A');
if min(eig(A)) > 0
    x0 = xc(1); y0 = xc(2);
    R = chol(A);
    theta = linspace(-pi,pi,100);
    xy_tilde = [cos(theta); sin(theta)];
    invR = inv(R);
    xy_bar = invR*xy_tilde;
    xe = xy_bar(1,:)+x0;
    ye = xy_bar(2,:)+y0;
    if strcmp(plotType,'patch')
        h = patch(xe,ye,color,'FaceAlpha',.2,'EdgeAlpha',0);
    else
        h = plot(xe,ye,'color',color, 'LineWidth',2);
    end
end
if min(eig(A))<=0
    %disp('invalid ellipsoid -- A not positive definite'), xe=[];, ye=[];, return
    lim = 8e4;
    x = linspace(-lim,lim,500);
    y = linspace(-lim,lim,500);
    a = A(1,1); b = A(1,2); c = A(2,2);
    yp = @(x) .5/b*(-2*c*x+sqrt(4*c^2*x^2 - 4*b*(a*x^2-1)));
    ym = @(x) .5/b*(-2*c*x-sqrt(4*c^2*x^2 - 4*b*(a*x^2-1)));
    plot(x, arrayfun(yp,x),[color linetype],'LineWidth',2);
    hold on
    plot(x, arrayfun(ym,x),[color linetype],'LineWidth',2);
    return
end
if ~type1
    xc=A\(-b./2);
    k=xc'*A*xc-c;
    if k<=0
        disp('invalid ellipsoid'), xe=[];, ye=[];, return
    end
    A=A/k;
end



