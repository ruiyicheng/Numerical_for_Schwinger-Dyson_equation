tic
high_temp=0.025;
low_temp=0.0005;
tstep=0.0005;
temperature=high_temp;
err=10^(-4);
global beta;
global lambda;
lambda=200;
beta=1/temperature;
j=1;
q=4;
x=0.5;
max=2*lambda;
mu=0.0;
omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
omega_n=gpuArray(omega_n);
g_tao=0*[-ones(1,2*lambda) ones(1,2*lambda)];
g_tao=gpuArray(g_tao);
temperature_list=[];
free_eng_list=[];
%----------------------------------------------------initialize
while temperature>=low_temp
    beta=1/temperature;
    dif=Inf;
    while dif>err*x
        g_omega=fourier(g_tao);
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
    temperature_list=[temperature_list temperature];
    p1=log(((1+s_omega./(1i.*omega_n)).^2+mu.^2./(omega_n).^2));
    p1=gather(sum(p1(lambda+1:3*lambda+1)));
    p2=real(gather(sum(s_omega(lambda+1:3*lambda+1).*g_omega(lambda+1:3*lambda+1))));
    p3=gather(sum(g_tao(2*lambda+1:4*lambda).^q));
    free_eng_list=[free_eng_list real(-log(2)./(2.* beta ) - p1/(4*beta)-1/(2.*beta)*p2-(p3*beta*j^2*2^(q-1))/(4*lambda*q^2))];
    temperature=temperature-tstep
    x=0.5;
end

temperature=low_temp

while temperature<high_temp
    beta=1/temperature;
    dif=Inf;
    while dif>err*x
        g_omega=fourier(g_tao);
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
    temperature_list=[temperature_list temperature];
    p1=log(((1+s_omega./(1i.*omega_n)).^2+mu.^2./(omega_n).^2));
    p1=gather(sum(p1(lambda+1:3*lambda+1)));
    p2=real(gather(sum(s_omega(lambda+1:3*lambda+1).*g_omega(lambda+1:3*lambda+1))));
    p3=gather(sum(g_tao(2*lambda+1:4*lambda).^q))
    free_eng_list=[free_eng_list real(-log(2)./(2.* beta ) - 1/(4*beta)*p1-1/(2.*beta)*p2-(beta*j^2*2^(q-1)*p3)/(4*lambda*q^2))];
    temperature=temperature+tstep
    x=0.5;
end
plot(temperature_list,2*free_eng_list)
tt=toc