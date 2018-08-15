% spherically uniform sample from a unit-sphere
% input: d=dimensionality
%        n=number of samples
% output:u=nxd matrix spherically uniform samples on a unit-sphere

function u=spherunifrnd(n,d)

mu=0;
sigma=1;
a=normrnd(mu,sigma,[n,d]);
b=sqrt(sum(a.^2,2));
b=repmat(b,1,d);
u=a./b;
   
% plot(u(:,1),u(:,2),'o')
end