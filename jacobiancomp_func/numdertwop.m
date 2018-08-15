% two-point method to compute numerical derivative
% inputs: fun=function handle of which derivative has to be computed
%         x=an array containing base points
% output: der=derivative computed by forward difference method

function der=numdertwop(fun,x)

der=(fun(x(2))-fun(x(1)))/(x(2)-x(1));

end