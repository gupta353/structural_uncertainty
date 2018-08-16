function [ pdf ] = GenGausspdf(param)
% This function computes pdf at input parameter
% for multivariate-Generalized Gaussian distribution
%%% inputs: 
% x=input column vector
%%% Output
% pdf=(log)probability density
%%% Reference: 
% Subbotin, (1923); On the law of frequency of errors.
% Gomez et al.,(1998) A multivariate generalization of power
% expoenntial family of distributions

global GLOBAL_DATA GEOMORPH

strmobs=GLOBAL_DATA.strmobs;

% parameters
theta=[param(1),param(2)];
a=[param(3),param(4)];
beta=param(5);
sigmaobs2=param(6);

% parameter perturbation covariance matrix
C=diag(a);

% computation of residuals
strmsim=hydrograph(theta,GLOBAL_DATA,GEOMORPH);
minlen=min(length(strmsim),length(strmobs));
strmsim=strmsim(1:minlen)';
strmobs_temp=strmobs(1:minlen);
err=errcompute(strmobs_temp,strmsim);  %% computation of error vector
% err(1)=[];
% 
n=length(err);
Ct=CovarianceMatrixEstimation(C,theta,sigmaobs2,GLOBAL_DATA,GEOMORPH);
log_fac=log(n)+log(gamma(n/2/beta))-(1/beta)*log(2)-log(gamma((n+2)/2/beta));
sigma=Ct;
sigma=sigma(1:n,1:n);
log_k=log(n)+log(gamma(n/2))-(n/2)*log(pi)-log(gamma(1+(n/2/beta)))-(1+(n/2/beta))*log(2);
sigmadet=det(sigma);

pdf=log_k-0.5*log(sigmadet)-0.5*n*log_fac-0.5*(err'*(sigma\err))^beta;
end