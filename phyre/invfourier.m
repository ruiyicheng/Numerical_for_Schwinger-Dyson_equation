function out = invfourier(inp)
% 1912.12302.B7
global lambda;
global beta;
lb=lambda;
%n=gpuArray(-lb:lb-1);
n=-lb:lb-1;
out=zeros(1,4*lambda);
%out=zeros(1,4*lambda,'gpuArray');
for m=-2*lb:2*lb-1
    idx=m+2*lb+1;
    aa=sum(exp(-(pi.*1i.*m.*(n + 0.5))./lb).*inp(lb+1:3*lb));
    out(idx)=aa;
end
out=out./beta;
end
