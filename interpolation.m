function [s, vals, valsf, iter] = interpolation(f, grad_f, x, d, s0)
%% title : Algorithme of quadratic interpolation combined to cubic interpolation
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% f       : Function to evaluate
% grad_f  : Function's gradient
% x       : Starting point
% d       : Descent direction
% s0      : Initial step
%% Output
% s       : Solution's approximation
% iter    : Number of iterations
% vals    : Step value
% valsf   : Values of f applied to x+s*d
%% 
global nfev ngev nhev;
c1 = 10^-3;
vals =[s0];
valsf=[f(x+s0*d)];
iter = 1;
if f(x+s0*d) <= f(x) + c1*s0*d'*grad_f(x)
    s=s0;
else
    s1 = -(d'*grad_f(x))/(2*((f(x+s0*d)-f(x)-s0*d'*grad_f(x))/s0^2));
    if (f(x+s1*d) <= f(x)+s1*c1*d'*grad_f(x))
        s=s1;
    else
        while ((f(x+s1*d) > f(x)+s1*c1*d'*grad_f(x)) && (abs(s1-s0)>1E-10) && iter<100) || iter==1
            iter = iter+1;
            cst =(s0^2*s1^2*(s1-s0));
            atemp = s0^2*(f(x+s1*d)-f(x)-d'*grad_f(x)*s1)-s1^2*(f(x+s0*d)-f(x)-d'*grad_f(x)*s0);
            btemp = -s0^3*(f(x+s1*d)-f(x)-d'*grad_f(x)*s1)+s1^3*(f(x+s0*d)-f(x)-d'*grad_f(x)*s0);
            a=(1/cst)*atemp;
            b=(1/cst)*btemp;
            delta = 4*b^2-12*a*d'*grad_f(x);
            r1 = (-2*b-sqrt(delta))/(6*a);
            r2 = (-2*b+sqrt(delta))/(6*a);
            s2 = max(r1,r2);
            if f(x+s2*d)<=f(x)+c1*s2*d'*grad_f(x)
                s=s2;
                break;
            else
                s0=s1;
                s1=s2;
            end
        end
        s=s2;
    end
end