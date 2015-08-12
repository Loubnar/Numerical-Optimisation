function [s, vals, valsf, iter] = approach(f, grad_f, x, d, s1, smax)
%% title : Approach Algorithm
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
% vals    : Values of the step
% valsf   : Values of f applied to x+s*d
%% 
global nfev ngev nhev;
c1=10^-3;
c2=0.9;
si=s1;
si_old=0;
s = 1;
vals = [s];
valsf = [f(x+s*d)];
iter = 0;
iterfinition = 0;
while iter<10000
    iter = iter+1;
    if (f(x+si*d) > f(x) + c1*si*d'*grad_f(x)) || (f(x+si*d)>=f(x+si_old*d) && iter>1)
        [s, iterfinition] = finition(f,grad_f,x,d,si_old,si,c1,c2);
        iter = iter + iterfinition;
        vals = [vals s];
        valsf = [valsf f(x+s*d)];
        break;
    elseif abs(d'*grad_f(x+si*d)) <= -c2*d'*grad_f(x)
        s = si;
        vals = [vals s];
        valsf = [valsf f(x+s*d)];
        break;
    elseif d'*grad_f(x+si*d)>=0
        s = finition(f,grad_f,x,d,si,si_old,c1,c2);
        vals = [vals s];
        valsf = [valsf f(x+s*d)];
        break;
    else
        si_old = si;
        si = si+(smax-si)*rand;
    end
end