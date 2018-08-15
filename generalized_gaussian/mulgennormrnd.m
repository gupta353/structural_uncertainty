% random sample draw from multivariate-generalized normal distribution
% inputs: n=number of samples
%         mu=mean (column vector)
%         cov=covraince matrix
%         beta=shape parameter
% output: r=nxd matrix with row cntaining a random sample, where d is the
% dimensionality of sample
% reference: Gomez et al., (1998). A multivariate generalization of the
% power exponential family of distributions

function r=mulgennormrnd(n,mu,cov,beta)

d=length(mu);

A=d*gamma(d/2/beta)/(2^(1/beta)*gamma((d+2)/2/beta))*chol(cov);

% spherically uniform sample
u=spherunifrnd(n,d);

% generation of radius random-variable
R=(gamrnd(d/2/beta,2,[1,n])).^(1/2/beta);      % check if it should be 2 or 1/2 (second argument in gamrnd)

r=repmat(mu,1,n)+(repmat(R,d,1)).*(A'*u');
r=r';

% plot(r(:,1),r(:,2),'o')
end