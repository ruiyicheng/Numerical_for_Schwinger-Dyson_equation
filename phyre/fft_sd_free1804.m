clear;
tic
%p=parpool;

global beta;
global lambda;

high_temp=0.0005;
low_temp=0.0005;
testep=0.0005;
temperature=high_temp;
err=10^(-15);

lambda=2^16; 
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



grlf_tao = -0.5*sign(points+0.5-(1:2*points));
gll_omega = -conj(0.5*tstep*fft(gll_tao));
glr_omega = -conj(0.5*tstep*fft(glr_tao));
grlf_omega = -conj(0.5*tstep*fft(gll_tao));
sll_tao = zeros(1,2*points);
sll_omega = zeros(1,2*points);
temperature_list=[];
free_eng_list=[];
wstep = pi/points;
w1 = zeros(1,2*points-1);
w1(1:points)=wstep*(1:points);
w1(points:2*points-1) = wstep*((points:2*points-1)-2*points);
delta=zeros(1,2*points);

%----------------------------------------------------initialize

while temperature>=low_temp
    dif=Inf;
    beta=1/temperature;
    delta(1)=-points/beta;
    delta(1+points)=points/beta;
    wstep = pi/points;
    w1 = zeros(1,2*points-1);
    w1(1:points)=wstep*(1:points);
    w1(points:2*points-1) = wstep*((points:2*points-1)-2*points);
    omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
    iter=0;
    tstep = beta/(2*lambda);
    while dif>err&&iter<5000
        %g_omega=parfourier(g_tao);
        sll_tao=j.*j.*(2.*gll_tao).^(q-1)./q;
        slr_tao=(-1)^(q/2)*j.*j.*(2.*glr_tao).^(q-1)./q-1i*mu.*delta;
        
        sll_omega=fftfourier(sll_tao);
        slr_omega=fftfourier(slr_tao);
        gll_omega(2:2*points)=((-(exp(-1i*w1)-1)/tstep+sll_omega(2:2*points))./((-(exp(-1i*w1)-1)/tstep+sll_omega(2:2*points)).^2+slr_omega(2:2*points).*fliplr(slr_omega(2:2*points)))).*mod(1:4*lambda-1,2);
        gll_tao_new=(1-x).*gll_tao+x.*fftinvfourier(gll_omega);
        
        glr_omega(2:2*points)=((-(exp(-1i*w1)-1)/tstep+sll_omega(2:2*points))./((-(exp(-1i*w1)-1)/tstep+sll_omega(2:2*points)).^2-mu ^2)).*mod(1:4*lambda-1,2);
        gll_tao_new=(1-x).*gll_tao+x.*fftinvfourier(gll_omega);
        nd=gather(sum(abs(gll_tao_new-gll_tao).^2));

        if nd<dif
            dif=nd;
            gll_tao=real(gll_tao_new);
            %plot(g_tao);
            %pause(0.01);
            iter=iter+1;
        else
            x=x/2;
            continue 
        end
    end
    grlf_omega(2:2*points)=((-1i*mu)./(((-(exp(-1i*w1)-1)/tstep+sll_omega(2:2*points)).^2-mu ^2))).*mod(1:4*lambda-1,2);
    grlf_tao=fftinvfourier(grlf_omega);
    grlf_omega=parfourier(grlf_tao);
    srlf_tao=(-1)^(q/2)*(j^2/q)*(2*grlf_tao).^(q-1)+1i*mu.*delta;
    srlf_omega=parfourier(srlf_tao);
    
    sllf_omega=parfourier(sll_tao);
    gllf_omega=(-(1i.*omega_n+sllf_omega)./((1i.*omega_n+sllf_omega).^2-mu^2));
    gllf_tao=gll_tao;
    sllf_tao=parinvfourier(sllf_omega);
    
    srrf_omega=sllf_omega;
    grrf_omega=gllf_omega;
    grrf_tao=gllf_tao;
    srrf_tao=sllf_tao;
    
    glrf_omega=-grlf_omega;
    glrf_tao=-grlf_tao;
    slrf_tao=(-1)^(q/2)*(j^2/q)*(2*glrf_tao).^(q-1)-1i*mu.*delta;
    slrf_omega=parfourier(slrf_tao);
    temperature_list=[temperature_list temperature];
    
    
    
    p1=2*log(((1+sllf_omega./(1i.*omega_n)).^2-(srlf_omega.*fliplr(srlf_omega))./(omega_n).^2));
    p1=gather(sum(p1(lambda+1:3*lambda+1)));
    p2=-beta^2/(4*lambda)*real(gather(sum(sllf_tao(2*lambda+1:4*lambda).*gllf_tao(2*lambda+1:4*lambda)+srrf_tao(2*lambda+1:4*lambda).*grrf_tao(2*lambda+1:4*lambda)+(slrf_tao(2*lambda+1:4*lambda).*glrf_tao(2*lambda+1:4*lambda)+srlf_tao(2*lambda+1:4*lambda).*grlf_tao(2*lambda+1:4*lambda)))));
    p3=gather(sum(gllf_tao(2*lambda+1:4*lambda).^q+(-1)^(q/2)*glrf_tao(2*lambda+1:4*lambda).^q+(-1)^(q/2)*grlf_tao(2*lambda+1:4*lambda).^q+grrf_tao(2*lambda+1:4*lambda).^q));
    p4=-1i*mu/2*(sum(grlf_omega(lambda+1:3*lambda+1)+glrf_omega(lambda+1:3*lambda+1)));
    free_eng_list=[free_eng_list real(-log(2)./(beta)-1/(4*beta)*p1-1/(beta)*p2-(beta*j^2*2^(q-1)*p3)/(4*lambda*q^2)+1/beta*p4)];
    temperature=temperature-testep
    x=0.5;
end
% 
temperature=low_temp
% 

while temperature<high_temp
    dif=Inf;
    beta=1/temperature;
    delta(1)=-points/beta;
    delta(1+points)=points/beta;
    wstep = pi/points;
    w1 = zeros(1,2*points-1);
    w1(1:points)=wstep*(1:points);
    w1(points:2*points-1) = wstep*((points:2*points-1)-2*points);
    omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
    iter=0;
    tstep = beta/(2*lambda);
    while dif>err&&iter<5000
        %g_omega=parfourier(g_tao);
        sll_tao=j.*j.*(2.*gll_tao).^(q-1)./q;
        sll_omega=fftfourier(sll_tao);
        gll_omega(2:2*points)=((-(exp(-1i*w1)-1)/tstep+sll_omega(2:2*points))./((-(exp(-1i*w1)-1)/tstep+sll_omega(2:2*points)).^2-mu ^2)).*mod(1:4*lambda-1,2);
        gll_tao_new=(1-x).*gll_tao+x.*fftinvfourier(gll_omega);
        nd=gather(sum(abs(gll_tao_new-gll_tao).^2));

        if nd<dif
            dif=nd;
            gll_tao=real(gll_tao_new);
            %plot(g_tao);
            %pause(0.01);
            iter=iter+1;
        else
            x=x/2;
            continue 
        end
    end
    grlf_omega(2:2*points)=((-1i*mu)./(((-(exp(-1i*w1)-1)/tstep+sll_omega(2:2*points)).^2-mu ^2))).*mod(1:4*lambda-1,2);
    grlf_tao=fftinvfourier(grlf_omega);
    grlf_omega=parfourier(grlf_tao);
    srlf_tao=(-1)^(q/2)*(j^2/q)*(2*grlf_tao).^(q-1)+1i*mu.*delta;
    srlf_omega=parfourier(srlf_tao);
    
    sllf_omega=parfourier(sll_tao);
    gllf_omega=(-(1i.*omega_n+sllf_omega)./((1i.*omega_n+sllf_omega).^2-mu^2));
    gllf_tao=gll_tao;
    sllf_tao=parinvfourier(sllf_omega);
    
    srrf_omega=sllf_omega;
    grrf_omega=gllf_omega;
    grrf_tao=gllf_tao;
    srrf_tao=sllf_tao;
    
    glrf_omega=-grlf_omega;
    glrf_tao=-grlf_tao;
    slrf_tao=(-1)^(q/2)*(j^2/q)*(2*glrf_tao).^(q-1)-1i*mu.*delta;
    slrf_omega=parfourier(slrf_tao);
    temperature_list=[temperature_list temperature];
    
    
    
    p1=2*log(((1+sllf_omega./(1i.*omega_n)).^2-(srlf_omega.*fliplr(srlf_omega))./(omega_n).^2));
    p1=gather(sum(p1(lambda+1:3*lambda+1)));
    p2=-beta^2/(4*lambda)*real(gather(sum(sllf_tao(2*lambda+1:4*lambda).*gllf_tao(2*lambda+1:4*lambda)+srrf_tao(2*lambda+1:4*lambda).*grrf_tao(2*lambda+1:4*lambda)+(slrf_tao(2*lambda+1:4*lambda).*glrf_tao(2*lambda+1:4*lambda)+srlf_tao(2*lambda+1:4*lambda).*grlf_tao(2*lambda+1:4*lambda)))));
    p3=gather(sum(gllf_tao(2*lambda+1:4*lambda).^q+(-1)^(q/2)*glrf_tao(2*lambda+1:4*lambda).^q+(-1)^(q/2)*grlf_tao(2*lambda+1:4*lambda).^q+grrf_tao(2*lambda+1:4*lambda).^q));
    p4=-1i*mu/2*(sum(grlf_omega(lambda+1:3*lambda+1)+glrf_omega(lambda+1:3*lambda+1)));
    free_eng_list=[free_eng_list real(-log(2)./(beta)-1/(4*beta)*p1-1/(beta)*p2-(beta*j^2*2^(q-1)*p3)/(4*lambda*q^2)+1/beta*p4)];    temperature=temperature+testep
    x=0.5;
end
hold off
plot(gll_tao);
hold on 
plot(imag(grlf_tao));
legend('g_tao','imag(glr_tao)')
tt=toc