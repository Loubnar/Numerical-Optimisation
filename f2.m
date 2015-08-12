function result = f2(x)
%% title : function f2
%% auteur : Loubna Rizqi
%% date : 01/2015
%% Input
% x       : Point
%% Output 
% result  : f1 applied to x
%%
global nfev;
nfev= nfev+1;
%% result
result = 100*(x(2)-x(1)^2)^2+(1-x(1))^2;
end