function out= parfourier(inp)
% 1912.12302.B7
global lambda;
global beta;
lb=lambda;
%m=gpuArray(0:2*lb-1);
m=0:2*lb-1;
%out=zeros(1,4*lambda,'gpuArray');
out=zeros(1,4*lambda);
ip=inp(2*lb+1:4*lb);
parfor k=1:4*lb
    aa=sum(exp((pi.*1i.*m.*(k -2*lb- 0.5))./lb).*ip);
    out(k)=aa;
end
out=out.*beta./(2*lambda);
end

