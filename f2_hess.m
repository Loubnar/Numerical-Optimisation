function result = f2_hess(x)
%% title : Hessian of function f2
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% x       : Point
%% Output 
% result  : f2 Hessian applied to x
%%
global nhev;
nhev = nhev+1;
%% result
result(1,:) = [-400*(x(2)-3*x(1)^2)+2 -400*x(1)];
result(2,:) = [-400*x(1) 200];
end

