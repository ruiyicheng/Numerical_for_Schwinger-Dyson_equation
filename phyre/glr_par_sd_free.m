clear;
tic
%p=parpool;

global beta;
global lambda;
high_temp=0.00086603;
low_temp=0.00086603;
tstep=0.0005;
temperature=high_temp;
err=10^(-1);

lambda=5000; 
beta=1/temperature;
j=1;
q=4;
x=0.3;
max=2*lambda;
mu=0.075;
omega_n= (2.*pi.*((1:2*max)-2*lambda- 0.5))./beta;
omega_n=gpuArray(omega_n);
gll_tao=0*[-ones(1,2*lambda) ones(1,2*lambda)];
gll_tao=gpuArray(gll_tao);
gll_omega=parfourier(gll_tao);
glr_tao=0*[-ones(1,2*lambda) ones(1,2*lambda)];
glr_tao=gpuArray(glr_tao);
glr_omega=parfourier(glr_tao);
delta=zeros(1,4*lambda);
delta(2*lambda+1)=2*lambda/beta;
delta(1)=-2*lambda/beta;
delta=gpuArray(delta);
temperature_list=[];
free_eng_list=[];
%----------------------------------------------------initialize




while temperature>=high_temp
    dif=Inf;
    beta=1/temperature;
    omega_n= (2.*pi.*((1:2*max)-max - 0.5))./beta;
    omega_n=gpuArray(omega_n);
    iter=0;
    while dif>err&&iter<150
        %g_omega=parfourier(g_tao);
        delta=zeros(1,4*lambda);
        delta(2*lambda+1)=-2*lambda/beta;
        delta(1)=2*lambda/beta;
        delta=gpuArray(delta);
        sll_tao=j.*j.*(2.*gll_tao).^(q-1)./q;
        slr_tao=j.*j.*(2.*glr_tao).^(q-1)./q.*(-1).^(q./2)-1i.*mu.*delta;
        
        sll_omega=parfourier(sll_tao);
        slr_omega=parfourier(slr_tao);
        
        gll_omega_new=(1-x).*gll_omega+x.*(1./(-1i.*omega_n-sll_omega+(slr_omega.^2/(-1i.*omega_n-sll_omega))));
        glr_omega_new=(1-x).*glr_omega+x.*((slr_omega.*gll_omega)./(-1i.*omega_n-sll_omega));
        
        gll_tao_new=parinvfourier(gll_omega_new);        
        glr_tao_new=parinvfourier(glr_omega_new);
        
        nd=gather(sum(abs(gll_omega_new-gll_omega).^2+abs(glr_omega_new-glr_omega).^2))

        if nd<gather(dif)
            dif=nd;
            gll_omega=gll_omega_new;
            gll_tao=real(gll_tao_new);
            glr_omega=glr_omega_new;
            glr_tao=1i.*imag(glr_tao_new);
            plot(gll_tao);
            hold on;
            plot(imag(glr_tao));
            hold off;
            pause(0.01)
            iter=iter+1
        else
            x=x/2;
            break;
            continue
        end
    end
    temperature_list=[temperature_list temperature];
    p1=log(((1+sll_omega./(1i.*omega_n)).^2+mu.^2./(omega_n).^2));
    p1=gather(sum(p1(lambda+1:3*lambda+1)));
    p2=real(gather(sum(sll_omega(lambda+1:3*lambda+1).*gll_omega(lambda+1:3*lambda+1))));
    p3=gather(sum(gll_tao(2*lambda+1:4*lambda).^q))
    free_eng_list=[free_eng_list real(-log(2)./(2.*beta)-1/(4*beta)*p1-1/(2.*beta)*p2-(beta*j^2*2^(q-1)*p3)/(4*lambda*q^2))];
    temperature=temperature-tstep
    x=0.5;
end
temperature=low_temp
while temperature<low_temp
    dif=Inf;
    beta=1/temperature;
    omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
    omega_n=gpuArray(omega_n);
    iter=0;
    while dif>err&&iter<150
        %g_omega=parfourier(g_tao);
        sll_tao=j.*j.*(2.*gll_tao).^(q-1)./q;
        sll_omega=parfourier(sll_tao);
        gll_omega_new=(1-x).*gll_omega-x.*((1i.*omega_n+sll_omega)./((1i.*omega_n+sll_omega).^2-mu^2));
        gll_tao_new=parinvfourier(gll_omega_new);
        nd=gather(sum(abs(gll_omega_new-gll_omega).^2))

        if nd<dif
            dif=nd;
            gll_omega=gll_omega_new;
            gll_tao=real(gll_tao_new);
            plot(gll_tao);
            pause(0.01);
            iter=iter+1
        else
            x=x/2;
            continue
        end
    end
    
    temperature_list=[temperature_list temperature];
    p1=log(((1+sll_omega./(1i.*omega_n)).^2+mu.^2./(omega_n).^2));
    p1=real(gather(sum(p1(lambda+1:3*lambda+1))));
    %p2=(gather(sum((s_omega(lambda+1:3*lambda+1)).*(g_omega(lambda+1:3*lambda+1)))));
    p2=-beta.^2./max.*(gather(sum((sll_tao(2*lambda+1:4*lambda)).*(gll_tao(2*lambda+1:4*lambda)))));
    p3=gather(sum((gll_tao(2*lambda+1:4*lambda)).^q));
    nf=(-log(2)./(2.*beta)-p1/(4*beta)-p2/(2.*beta)-(beta*(j^2)*2.^(q-1))/(4*lambda*q^2)*p3);
    free_eng_list=[free_eng_list nf];
    temperature=temperature+tstep
    x=0.5;
end




scatter(temperature_list,2*free_eng_list)
tt=toc