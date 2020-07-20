clear;
tic
%p=parpool;

global beta;
global lambda;

high_temp=0.04;
low_temp=0.001;
testep=0.001;
temperature=high_temp;
err=10^(-5);

lambda=2^8; 
beta=1/temperature;
j=1;
q=4;
x=0.5;
max=2*lambda;
points=max;
mu=0.03;
omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
tstep = beta/(2*lambda);

gll_tao = -0.5*sign(points+0.5-(1:2*points));
glr_tao=zeros(1,2*points);
grr_tao = -0.5*sign(points+0.5-(1:2*points));
grl_tao=zeros(1,2*points);




gll_omega = parfourier(gll_tao);
glr_omega =parfourier(glr_tao);
grr_omega = parfourier(grr_tao);
grl_omega =parfourier(grl_tao);


temperature_list=[];
free_eng_list=[];

delta=zeros(1,2*points);

%----------------------------------------------------initialize

while temperature>=low_temp
    dif=Inf;
    beta=1/temperature;
    delta(1)=-points/beta;
    delta(1+points)=points/beta;

    omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
    iter=0;
    tstep = beta/(2*lambda);
    while dif>err&&iter<5000
        %g_omega=parfourier(g_tao);
        sll_tao=j.*j.*(2.*gll_tao).^(q-1)./q;
        slr_tao=(-1)^(q/2)*j.*j.*(2.*glr_tao).^(q-1)./q-1i*mu.*delta;
        srr_tao=j.*j.*(2.*grr_tao).^(q-1)./q;
        srl_tao=(-1)^(q/2)*j.*j.*(2.*grl_tao).^(q-1)./q+1i*mu.*delta;
        
        
        sll_omega=parfourier(sll_tao);
        slr_omega=parfourier(slr_tao);
        srr_omega=parfourier(srr_tao);
        srl_omega=parfourier(srl_tao);        
        
        A=-slr_omega.*srl_omega+(-(exp(1i*omega_n*tstep)-1)./tstep-srr_omega).*(-(exp(1i*omega_n*tstep)-1)./tstep-sll_omega);
       
        
        gll_omega=(-(exp(1i*omega_n*tstep)-1)./tstep-srr_omega)./A;
        glr_omega=slr_omega./A;
        grl_omega=srl_omega./A;
        grr_omega=(-(exp(1i*omega_n*tstep)-1)./tstep-sll_omega)./A;

        
        
        gll_tao_new=(1-x).*gll_tao+x.*parinvfourier(gll_omega);
        glr_tao_new=(1-x).*glr_tao+x.*parinvfourier(glr_omega);    
        grr_tao_new=(1-x).*grr_tao+x.*parinvfourier(grr_omega);        
        grl_tao_new=(1-x).*grl_tao+x.*parinvfourier(grl_omega);
        
        nd=gather(sum(abs(gll_tao_new-gll_tao).^2))

        if nd<dif
            dif=nd;
%             gll_omega=gll_omega_new;
%             glr_omega=glr_omega_new;
%             grl_omega=grl_omega_new;
%             grr_omega=grr_omega_new;
            gll_tao=real(gll_tao_new);
            glr_tao=1i*imag(glr_tao_new);
            grr_tao=real(grr_tao_new);
            grl_tao=1i*imag(grl_tao_new);
            hold off
            plot(gll_tao);
            hold on 
            plot(imag(glr_tao));
            plot(real(grr_tao));
            plot(imag(grl_tao));
            pause(0.01);
            iter=iter+1;
        else
            x=x/2
            continue 
        end
    end
     
     slrf_tao=slr_tao;
     glrf_tao=glr_tao;
     slrf_omega=parfourier(slr_tao);
     glrf_omega=parfourier(glr_tao);
     
%     
     sllf_tao=sll_tao;
     gllf_tao=gll_tao;
     sllf_omega=parfourier(sll_tao);
     gllf_omega=parfourier(gll_tao);
%     
     srrf_tao=srr_tao;
     grrf_tao=grr_tao;
     srrf_omega=parfourier(srr_tao);
     grrf_omega=parfourier(grr_tao);
%     
     srlf_tao=srl_tao;
     grlf_tao=grl_tao;
     srlf_omega=parfourier(srl_tao);
     grlf_omega=parfourier(grl_tao);
     
     temperature_list=[temperature_list temperature];
%     
%     
%     
     Seff1=-0.5*log(((1+sllf_omega./(1i.*omega_n)).*(1+srrf_omega./(1i.*omega_n))+(slrf_omega.*srlf_omega)./(omega_n).^2));
     Seff1=real(sum(Seff1(lambda+1:3*lambda+1))-log(2));
     
     Seff2=beta/(4*lambda)*real(sum(sllf_tao(2*lambda+1:4*lambda).*gllf_tao(2*lambda+1:4*lambda)+srrf_tao(2*lambda+1:4*lambda).*grrf_tao(2*lambda+1:4*lambda)+slrf_tao(2*lambda+1:4*lambda).*glrf_tao(2*lambda+1:4*lambda)+srlf_tao(2*lambda+1:4*lambda).*grlf_tao(2*lambda+1:4*lambda)));
     
     Seff3=-beta^2*j^2*2^q/(8*lambda*q^2)*real(sum(gllf_tao(2*lambda+1:4*lambda).^q+(-1)^(q/2)*glrf_tao(2*lambda+1:4*lambda).^q+(-1)^(q/2)*grlf_tao(2*lambda+1:4*lambda).^q+grrf_tao(2*lambda+1:4*lambda).^q));
 
     Seff4=real(1i*mu/2*sum(-grlf_omega(lambda+1:3*lambda+1)+glrf_omega(lambda+1:3*lambda+1)));
     
     free_eng_list=[free_eng_list real(Seff1+Seff2+Seff3+Seff4)/beta];
    temperature=temperature-testep
    x=0.5;
end
% 
temperature=low_temp
% 
% 
% omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
% tstep = beta/(2*lambda);
% gll_tao = -0.5*sign(points+0.5-(1:2*points));
% glr_tao=zeros(1,2*points);
% 
% 
% 
% grlf_tao = -0.5*sign(points+0.5-(1:2*points));
% gll_omega = -conj(0.5*tstep*fft(gll_tao));
% glr_omega = -conj(0.5*tstep*fft(glr_tao));
% grlf_omega = -conj(0.5*tstep*fft(gll_tao));
% sll_tao = zeros(1,2*points);
% slr_tao = zeros(1,2*points);
% sll_omega = zeros(1,2*points);
% 
% wstep = pi/points;
% w1 = zeros(1,2*points-1);
% w1(1:points)=wstep*(1:points);
% w1(points:2*points-1) = wstep*((points:2*points-1)-2*points);
% delta=zeros(1,2*points);

while temperature<high_temp
    dif=Inf;
    beta=1/temperature;
    delta(1)=-points/beta;
    delta(1+points)=points/beta;

    omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
    iter=0;
    tstep = beta/(2*lambda);
    while dif>err&&iter<5000
        %g_omega=parfourier(g_tao);
        sll_tao=j.*j.*(2.*gll_tao).^(q-1)./q;
        slr_tao=(-1)^(q/2)*j.*j.*(2.*glr_tao).^(q-1)./q-1i*mu.*delta;
        srr_tao=j.*j.*(2.*grr_tao).^(q-1)./q;
        srl_tao=(-1)^(q/2)*j.*j.*(2.*grl_tao).^(q-1)./q+1i*mu.*delta;
        
        
        sll_omega=parfourier(sll_tao);
        slr_omega=parfourier(slr_tao);
        srr_omega=parfourier(srr_tao);
        srl_omega=parfourier(srl_tao);        
        
        A=-slr_omega.*srl_omega+(-(exp(1i*omega_n*tstep)-1)./tstep-srr_omega).*(-(exp(1i*omega_n*tstep)-1)./tstep-sll_omega);
       
        
        gll_omega=(-(exp(1i*omega_n*tstep)-1)./tstep-srr_omega)./A;
        glr_omega=slr_omega./A;
        grl_omega=srl_omega./A;
        grr_omega=(-(exp(1i*omega_n*tstep)-1)./tstep-sll_omega)./A;

        
        
        gll_tao_new=(1-x).*gll_tao+x.*parinvfourier(gll_omega);
        glr_tao_new=(1-x).*glr_tao+x.*parinvfourier(glr_omega);    
        grr_tao_new=(1-x).*grr_tao+x.*parinvfourier(grr_omega);        
        grl_tao_new=(1-x).*grl_tao+x.*parinvfourier(grl_omega);
        
        nd=gather(sum(abs(gll_tao_new-gll_tao).^2))

        if nd<dif
            dif=nd;
%             gll_omega=gll_omega_new;
%             glr_omega=glr_omega_new;
%             grl_omega=grl_omega_new;
%             grr_omega=grr_omega_new;
            gll_tao=real(gll_tao_new);
            glr_tao=1i*imag(glr_tao_new);
            grr_tao=real(grr_tao_new);
            grl_tao=1i*imag(grl_tao_new);
            hold off
            plot(gll_tao);
            hold on 
            plot(imag(glr_tao));
            plot(real(grr_tao));
            plot(imag(grl_tao));
            pause(0.01);
            iter=iter+1;
        else
            x=x/2
            continue 
        end
    end
     
     slrf_tao=slr_tao;
     glrf_tao=glr_tao;
     slrf_omega=parfourier(slr_tao);
     glrf_omega=parfourier(glr_tao);
     
%     
     sllf_tao=sll_tao;
     gllf_tao=gll_tao;
     sllf_omega=parfourier(sll_tao);
     gllf_omega=parfourier(gll_tao);
%     
     srrf_tao=srr_tao;
     grrf_tao=grr_tao;
     srrf_omega=parfourier(srr_tao);
     grrf_omega=parfourier(grr_tao);
%     
     srlf_tao=srl_tao;
     grlf_tao=grl_tao;
     srlf_omega=parfourier(srl_tao);
     grlf_omega=parfourier(grl_tao);
     
     temperature_list=[temperature_list temperature];
%     
%     
%     
     Seff1=-0.5*log(((1+sllf_omega./(1i.*omega_n)).*(1+srrf_omega./(1i.*omega_n))+(slrf_omega.*srlf_omega)./(omega_n).^2));
     Seff1=real(sum(Seff1(lambda+1:3*lambda+1))-log(2));
     
     Seff2=beta/(4*lambda)*real(sum(sllf_tao(2*lambda+1:4*lambda).*gllf_tao(2*lambda+1:4*lambda)+srrf_tao(2*lambda+1:4*lambda).*grrf_tao(2*lambda+1:4*lambda)+slrf_tao(2*lambda+1:4*lambda).*glrf_tao(2*lambda+1:4*lambda)+srlf_tao(2*lambda+1:4*lambda).*grlf_tao(2*lambda+1:4*lambda)));
     
     Seff3=-beta^2*j^2*2^q/(8*lambda*q^2)*real(sum(gllf_tao(2*lambda+1:4*lambda).^q+(-1)^(q/2)*glrf_tao(2*lambda+1:4*lambda).^q+(-1)^(q/2)*grlf_tao(2*lambda+1:4*lambda).^q+grrf_tao(2*lambda+1:4*lambda).^q));
 
     Seff4=real(1i*mu/2*sum(-grlf_omega(lambda+1:3*lambda+1)+glrf_omega(lambda+1:3*lambda+1)));
     
     free_eng_list=[free_eng_list real(Seff1+Seff2+Seff3+Seff4)/beta];
    temperature=temperature+testep
    x=0.5;
end
hold off
% plot(gll_tao);
% hold on 
% plot(imag(glr_tao));
% legend('g_tao','imag(glr_tao)')
plot(temperature_list, free_eng_list)
tt=toc