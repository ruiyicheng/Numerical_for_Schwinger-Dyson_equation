function out= fourier(inp)
% 1912.12302.B7
global lambda;
global beta;
lb=lambda;
m=gpuArray(0:2*lb-1);

out=zeros(1,4*lambda,'gpuArray');

for k=-2*lb:2*lb-1
    idx=k+2*lb+1;
    aa=sum(exp((pi.*1i.*m.*(k + 0.5))./lb).*inp(2*lb+1:4*lb));
    out(idx)=aa;
end
out=out.*beta./(2*lambda);
end

