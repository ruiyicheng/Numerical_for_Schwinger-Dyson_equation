function out= fftinvfourier(inp)
% 1912.12302.B7
global lambda;
global beta;
lb=lambda;
out=4*lb*ifft(conj(inp))/beta;
end

