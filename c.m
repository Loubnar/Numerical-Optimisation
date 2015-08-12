function ceq = c(x)
global nfcev;
nfcev = nfcev+1;
%Computing the constraints c

ceq = (x(1)^2+x(2)^2)-1.5;
end