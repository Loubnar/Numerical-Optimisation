function [s, vals, valsf, exiter] = bissection(f, grad_f, x, d)
%% title : Bissection Algorithm
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% f       : Function to evaluate
% grad_f  : Function's gradient
% x       : Starting point
% epsilon : Precision of criteria
%% Output 
% s       : Solution's approximation
% iter    : Number of iterations
% vals    : Step value
% valsf   : Values of f applied to x+s*d
%% 
global nfev ngev nhev;
beta = inf;
alpha = 0;
c1 = 10^-3;
c2 = 0.9;
s = 1;
vals = [s];
valsf = [f(x+s*d)];
exiter = 0;
while true
    exiter = exiter + 1;
    if f(x+s*d) > f(x) + s*c1*grad_f(x)'*d
        % If the first condition of Wolfe is not met...
        beta = s;
        s = 0.5*(alpha+beta);
    elseif grad_f(x+s*d)'*d < c2*grad_f(x)'*d
        % If the second condition of Wolfe is not met...
        alpha = s;
        if beta == inf
            s = 2*alpha;
        else
            s = 0.5*(alpha+beta);
        end
    else
        break;
    end
    vals = [vals s];
    valsf = [valsf f(x+s*d)];
end

