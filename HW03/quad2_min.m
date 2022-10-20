function [x,y] = quad2_min(eps,alpha,beta)
% [x,y] = QUAD2_MIN(eps,alpha,beta)
%
% Solves:
%      minimize    (x+y)^2+eps*(x^2+y^2)+alpha*x+beta*y
%      subject to  x>=0, y>=0

x = (-1/((2*eps*(2+eps))))*((1+eps)*alpha-beta);
y = (-1/((2*eps*(2+eps))))*((1+eps)*beta-alpha);
        
if(x<0)||(y<0)
    if alpha<beta
        x = max(-alpha/(2*(1+eps)),0);
        y = 0;
    elseif alpha>beta
        y = max(-beta/(2*(1+eps)),0);
        x = 0;
    elseif alpha==beta
        x = max(-alpha/(2*(2+eps)),0);
        y = max(-alpha/(2*(2+eps)),0);
    end
end