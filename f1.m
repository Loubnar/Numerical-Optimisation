function result = f1(x)
%% title : function f1
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% x       : Point
%% Output 
% result  : f1 applied to w
%%
global nfev;
nfev = nfev+1;
%% Result
result = 2*(x(1)+x(2)-2)^2+(x(1)-x(2))^2;
end

