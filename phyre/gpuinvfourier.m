function out = gpuinvfourier(inp)
% 1912.12302.B7
global lambda;
global beta;
b=beta;
lb=lambda;
n=gpuArray((-lb:lb-1)');
m=gpuArray((1:4*lb));
out=zeros(1,4*lambda,'gpuArray');
ip=inp(lb+1:3*lb);
out=ip*exp(-(pi.*1i.*(m-2*lb-1).*(n + 0.5))./lb);
%parfor m=1:4*lb
%    aa=
%    (m)=aa;
%end
out=out./b;
end
