function result = J_c( x )
global ngcev;
ngcev =ngcev+1;
% Computing the jacobian of the constraints c

result = [2*x(1);2*x(2)];

end