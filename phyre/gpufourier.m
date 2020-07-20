function out= gpufourier(inp)
% 1912.12302.B7
global lambda;
global beta;
b=beta;
lb=lambda;
m=gpuArray((0:2*lb-1)');
k=gpuArray(1:4*lb);
out=zeros(1,4*lambda,'gpuArray');
ip=inp(2*lb+1:4*lb);
%sum(exp((pi.*1i.*m.*(k-2*lb- 0.5))./lb));
out=ip*exp((pi.*1i.*m.*(k-2*lb- 0.5))./lb);
%parfor k=1:4*lb
%    aa=sum(exp((pi.*1i.*m.*(k-2*lb- 0.5))./lb).*ip);
%    out(k)=aa;
%end
out=out.*b./(2*lambda);
end

