function [s, iter] = finition(f, grad_f, x, d, smin, smax, c1, c2)
%% title : Finition Algorithm
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% f       : Function to evaluate
% grad_f  : Function's gradient
% x       : Starting point
% d       : Descent direction
% s1      : Initial step
% smax    : Maximum step
% c1      : coefficient of the first condition of Wolfe
% c2      : coefficient of the second condition of Wolfe
%% Output
% s       : Satisfying step
% iter    : Number of iterations
%% 
global nfev ngev nhev;
iter = 0;
while iter<50000
    iter =iter+1;
    sj = smin+(smax-smin)*rand;
    if (f(x+sj*d) > f(x)+c1*sj*d'*grad_f(x)) && f(x+sj*d)>=f(x+smin*d)
        smax = sj;
    else
        if abs(d'*grad_f(x+sj*d))<= -c2*d'*grad_f(x)
            s = sj;
            break;
        elseif d'*grad_f(x+sj*d)*(smax-smin)>=0
            smax = smin;
            smin = sj;
        else
            smin = sj;
        end
    end          
end
s=sj;  