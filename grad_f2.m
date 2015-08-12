function result = grad_f2(x)
%% title : Gradient of f2
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% x       : Point
%% Output 
% result  : Gradient of f2 applied to x
%%
global ngev;
ngev = ngev+1;
%% result
result = [-400*x(1)*(x(2)-x(1)^2)-2*(1-x(1));200*(x(2)-x(1)^2)];
end

