% objective function defined as negative-log likelihood of the
% observations; the elements of covariance matrix and the model parameters
% are to be estimated
% input: param: param(1)=vel_stream (velocity of water drops in streams)
%               param(2)=vel_hillslope (velocity of water drops in hillslope)
%               param(3)=variance of perturbation of vel_stream
%               param(4)=variance of perturbation of vel_hillslope
%               param(5)=shape(kurtosis) parameter
%               param(6)=variance of observation errors
% output: negative-log likelihood evaluated at given 'parameters' for
% multivariate-generalized Gaussian distribution
% note: the observations which were not sesitive to the parameters are
% removed in computation of objective function. For this purpose, only the
% obsevation with computed variance greater than 0.001 are kept.
% reference: Gomez et al. (1998). A multivariate generalization of the
% power exponential family of distributions
function L = objectiveComputation(param,GLOBAL_DATA,GEOMORPH)

strmobs=GLOBAL_DATA.strmobs;

% parameters to be optimized
theta=[param(1),param(2)];
a=[param(3),param(4)];
beta=param(5);
sigmaobs2=param(6);             % observation variance

% computation of residuals
strmsim=hydrograph(theta,GLOBAL_DATA,GEOMORPH);
minlen=min(length(strmsim),length(strmobs));
strmsim=strmsim(1:minlen)';
strmobs_temp=strmobs(1:minlen);
err=errcompute(strmobs_temp,strmsim);  %% computation of error vector
% err(1)=[];
% definition of parameter-perturbation covariance matrix
C = diag(a);

% computation of total-error covariance matrix
Ct=CovarianceMatrixEstimation(C,theta,sigmaobs2,GLOBAL_DATA,GEOMORPH);
n=length(err); % number of components in err-vector
log_fac=log(n)+log(gamma(n/2/beta))-(1/beta)*log(2)-log(gamma((n+2)/2/beta));
sigma=Ct;
log_k=log(n)+log(gamma(n/2))-(n/2)*log(pi)-log(gamma(1+(n/2/beta)))-(1+(n/2/beta))*log(2);
sigma=sigma(1:n,1:n);
sigmadet=det(sigma);
L=-log_k+0.5*log(sigmadet)+0.5*n*log_fac+0.5*(err'*(sigma\err))^(beta);

end