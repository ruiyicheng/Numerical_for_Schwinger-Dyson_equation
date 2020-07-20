function out= fftfourier(inp)
% 1912.12302.B7
global lambda;
global beta;
lb=lambda;
tstep = beta/(2*lb);
out=-conj(0.5*tstep*fft(inp));
end

