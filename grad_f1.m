function result = grad_f1(x)
%% title : Gradient of f1
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% x       : Point
%% Output 
% result  : Gradient of f1 applied to x
%%
global ngev;
ngev = ngev+1;
%% result
result = [4*(x(1)+x(2)-2)+2*(x(1)-x(2));4*(x(1)+x(2)-2)-2*(x(1)-x(2))];
end

