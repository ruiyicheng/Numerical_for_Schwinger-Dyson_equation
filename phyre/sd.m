tic
temperature=0.04;
err=10^(-2);
global beta;
global lambda;
lambda=250;
beta=1/temperature;
j=1;
q=4;
x=0.5;
max=2*lambda;
mu=0.075;
omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
omega_n=gpuArray(omega_n);
g_tao=0*[-ones(1,2*lambda) ones(1,2*lambda)];
g_tao=gpuArray(g_tao);
dif=Inf;%----------------------------------------------------initialize

while dif>err*x
    %g_omega=fourier(g_tao);
    s_tao=j.*j.*(2.*g_tao).^(q-1)./q;
    s_omega=fourier(s_tao);
    g_omega_new=(1-x).*g_omega+x.*(-(1i.*omega_n+s_omega)./((1i.*omega_n+s_omega).^2-mu^2));
    g_tao_new=invfourier(g_omega_new);
    nd=gather(sum(abs(g_tao_new-g_tao).^2))
    
    if nd<gather(dif)
        dif=nd;
        g_omega=g_omega_new;
        g_tao=real(g_tao_new);
    else
        x=x/2
        continue
    end
end
plot(real(gather(g_tao)));

tt=toc