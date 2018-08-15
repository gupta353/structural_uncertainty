% computation of the covariance matrix
% input: C=covariance matrix of parameter perturbation
%        x=hydrologic-model parameter
%        sigmaobs2=variance of observation error
% output: F=covariance matrix of total error

function [F] = CovarianceMatrixEstimation(C,theta,sigmaobs2,GLOBAL_DATA,GEOMORPH)

J = jacobiancomp(@(x)hydrograph(x,GLOBAL_DATA,GEOMORPH),theta,0.01)';
F=J*C*J'+sigmaobs2*eye(size(J,1));
    
end

