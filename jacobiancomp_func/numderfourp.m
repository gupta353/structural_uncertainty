% four point method for computing numerical derivatives
% inputs: fun=function handle of which derivative has to be computed
%         a=the point at which derivative has to be computed
%         h=step-size for selcting four points
% output: der=derivative of 'fun' at 'a'
% Reference: Numerical differentiation and integration (Chapter 11)

function der=numderfourp(fun,x)

der=(fun(x(1))-8*fun(x(2))+8*fun(x(3))-fun(x(4)))/(12*(x(2)-x(1)));

end

