function out = parinvfourier(inp)
% 1912.12302.B7
global lambda;
global beta;
lb=lambda;
%n=gpuArray(-lb:lb-1);
n=(-lb:lb-1);
%out=zeros(1,4*lambda,'gpuArray');
out=zeros(1,4*lambda);
ip=inp(lb+1:3*lb);
parfor m=1:4*lb
    aa=sum(exp(-(pi.*1i.*(m-1-2*lb).*(n + 0.5))./lb).*ip);
    out(m)=aa;
end
out=out./beta;
end
