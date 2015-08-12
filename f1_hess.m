function result = f1_hess(x)
%% title : Hessian of f1
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% x       : Point
%% Output 
% result  : f1 hessian applied to x
%%
global nhev;
nhev = nhev+1;
%% result
result(1,:) = [6 2];
result(2,:) = [2 6];
end

