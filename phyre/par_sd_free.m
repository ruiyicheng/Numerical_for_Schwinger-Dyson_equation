clear;
tic
%p=parpool;

global beta;
global lambda;
high_temp=0.000866;
low_temp=0.000866;
tstep=0.0005;
temperature=high_temp;
err=10^(-3);

lambda=50000; 
beta=1/temperature;
j=1;
q=4;
x=0.5;
max=2*lambda;
mu=0.75;
omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
omega_n=gpuArray(omega_n);
g_tao=0*[-ones(1,2*lambda) ones(1,2*lambda)];
g_tao=gpuArray(g_tao);
g_omega=parfourier(g_tao);
temperature_list=[];
free_eng_list=[];
%----------------------------------------------------initialize

while temperature>=low_temp
    dif=Inf;
    beta=1/temperature;
    omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
    omega_n=gpuArray(omega_n);
    iter=0;
    while dif>err&&iter<150
        %g_omega=parfourier(g_tao);
        s_tao=j.*j.*(2.*g_tao).^(q-1)./q;
        s_omega=parfourier(s_tao);
        g_omega_new=(1-x).*g_omega-x.*((1i.*omega_n+s_omega)./((1i.*omega_n+s_omega).^2-mu^2));
        g_tao_new=parinvfourier(g_omega_new);
        nd=gather(sum(abs(g_omega_new-g_omega).^2))

        if nd<dif
            dif=nd;
            g_omega=g_omega_new;
            g_tao=real(g_tao_new);
            plot(g_tao);
            pause(0.01);
            iter=iter+1
        else
            x=x/2;
            continue
        end
    end
    
    temperature_list=[temperature_list temperature];
    p1=log(((1+s_omega./(1i.*omega_n)).^2+mu.^2./(omega_n).^2));
    p1=real(gather(sum(p1(lambda+1:3*lambda+1))));
    %p2=(gather(sum((s_omega(lambda+1:3*lambda+1)).*(g_omega(lambda+1:3*lambda+1)))));
    p2=-beta.^2./max.*(gather(sum((s_tao(2*lambda+1:4*lambda)).*(g_tao(2*lambda+1:4*lambda)))));
    p3=gather(sum((g_tao(2*lambda+1:4*lambda)).^q));
    nf=(-log(2)./(2.*beta)-p1/(4*beta)-p2/(2.*beta)-(beta*(j^2)*2.^(q-1))/(4*lambda*q^2)*p3);
    free_eng_list=[free_eng_list nf];
    temperature=temperature-tstep
    x=0.5;
end

temperature=low_temp

while temperature<high_temp
    dif=Inf;
    beta=1/temperature;
    omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
    omega_n=gpuArray(omega_n);
    iter=0
    while dif>err&&iter<150
        %g_omega=parfourier(g_tao);
        s_tao=j.*j.*(2.*g_tao).^(q-1)./q;
        s_omega=parfourier(s_tao);
        g_omega_new=(1-x).*g_omega+x.*(-(1i.*omega_n+s_omega)./((1i.*omega_n+s_omega).^2-mu^2));
        g_tao_new=parinvfourier(g_omega_new);
        nd=gather(sum(abs(g_omega_new-g_omega).^2))

        if nd<gather(dif)
            dif=nd;
            g_omega=g_omega_new;
            g_tao=real(g_tao_new);
            plot(g_tao);
            pause(0.01)
            iter=iter+1
        else
            x=x/2;
            continue
        end
    end
    temperature_list=[temperature_list temperature];
    p1=log(((1+s_omega./(1i.*omega_n)).^2+mu.^2./(omega_n).^2));
    p1=gather(sum(p1(lambda+1:3*lambda+1)));
    p2=real(gather(sum(s_omega(lambda+1:3*lambda+1).*g_omega(lambda+1:3*lambda+1))));
    p3=gather(sum(g_tao(2*lambda+1:4*lambda).^q))
    free_eng_list=[free_eng_list real(-log(2)./(2.*beta)-1/(4*beta)*p1-1/(2.*beta)*p2-(beta*j^2*2^(q-1)*p3)/(4*lambda*q^2))];
    temperature=temperature+tstep
    x=0.5;
end





scatter(temperature_list,2*free_eng_list)
tt=toc