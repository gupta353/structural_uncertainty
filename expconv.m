% convolution of independent nonidentical exponential pdfs
% input: lambda = parameters (1/mean-time) of each pdf to be convolved
%                 (row vector) with each element of array conatining
%                 one parameter
% output: uh_fun = a cell array with a function handle in each cell
%         the convolved value at each point the sum of functions in each
%         cell

%%%%%%%%%%%%%%%%%%% the mathematical form of pdf assumed %%%%%%%%%%%%%%%%%%
% f(x)=lambda*exp(-lambda*x)*I_(0,inf)(x),
% where lambda= parameter (1/mean-time) and I=indicator function

function uh_fun=expconv(lambda)
    
    k=length(lambda);
%     Lambda=repmat(lambda,k,1);
    temp_den=bsxfun(@minus,lambda',lambda);
    ind=temp_den==0;
    
    if sum(ind(:))>k
        error('atleast two of the pdfs to be convolved are identical')
    else
        temp_den(ind)=1;
        den=prod(temp_den);
        nume=prod(lambda);
        C=nume./den;
        uh_fun=cell(1,k);
        for i=1:k
            uh_fun{i}=@(t)exp(-lambda(i)*t)*C(i);
        end
    end
    
end

