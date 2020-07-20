clear;
tic
%p=parpool;

global beta;
global lambda;

high_temp=0.04;
low_temp=0.001;
testep=0.001;
temperature=high_temp;
err=10^(-15);
mu=0.03;
kappa=0.04;
lambda=2^11; 
beta=1/temperature;
J=1;
x=0.5;
max=2*lambda;
points=max;
omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
tstep = beta/(2*lambda);
g11_tao = -0.5*sign(points+0.5-(1:2*points));
g12_tao=zeros(1,2*points);



grlf_tao = -0.5*sign(points+0.5-(1:2*points));
g11_omega = -conj(0.5*tstep*fft(g11_tao));
g12_omega = -conj(0.5*tstep*fft(g12_tao));
grlf_omega = -conj(0.5*tstep*fft(g11_tao));
s11_tao = zeros(1,2*points);
s11_omega = zeros(1,2*points);
temperature_list=[];
free_eng_list=[];
wstep = pi/points;
w1 = zeros(1,2*points-1);
w1(1:points)=wstep*(1:points);
w1(points:2*points-1) = wstep*((points:2*points-1)-2*points);


%----------------------------------------------------initialize

while temperature>=low_temp
    dif=Inf;
    beta=1/temperature;
    wstep = pi/points;
    w1 = zeros(1,2*points-1);
    w1(1:points)=wstep*(1:points);
    w1(points:2*points-1) = wstep*((points:2*points-1)-2*points);
    omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
    iter=0;
    tstep = beta/(2*lambda);
    while dif>err&&iter<5000
        %g_omega=parfourier(g_tao);
        s11_tao=-J.*J.*g11_tao.^2.*fliplr(g11_tao);
        s12_tao=J.*J.*g12_tao.^2.*fliplr(g12_tao);
        
        s11_omega=parfourier(s11_tao);
        s12_omega=parfourier(s12_tao);
        D=(-(exp(1i.*omega_n*tstep)-1)/tstep-mu-s11_omega).^2+(1i.*kappa-s12_omega).^2;
        
        g11_omega=(-(exp(1i.*omega_n*tstep)-1)/tstep-mu-s11_omega)./D;
        gll_tao_new=(1-x).*g11_tao+x.*parinvfourier(g11_omega);
        
        g12_omega=(-1i.*kappa+s12_omega)./D;
        gl2_tao_new=(1-x).*g12_tao+x.*parinvfourier(g12_omega);
        nd=gather(sum(abs(gll_tao_new-g11_tao).^2));

        if nd<dif
            dif=nd;
            g11_tao=real(gll_tao_new);
            g12_tao=1i*imag(gl2_tao_new);
            plot(g11_tao);
            hold on
            plot(-imag(g12_tao))
            hold off
            pause(0.1);
            iter=iter+1;
        else
            x=x/2;
            continue 
        end
    end
     D=(-1i.*omega_n-mu-s11_omega).^2+(1i.*kappa-s12_omega).^2;
%     D=(-1i.*omega_n-mu-s11_omega).^2+(1i.*kappa-s12_omega).^2;
%     s11f_omega=parfourier(s11_tao);
%     g11f_omega=(-1i.*omega_n-mu-s11f_omega)./D;
%     g11f_tao=g11_tao;
%     s11f_tao=parinvfourier(s11f_omega);
%     
% 
%     g12f_omega=parfourier(s12_tao);
%     g12f_tao=g12_tao;
%     s12f_tao=s12_tao;
%     s12f_omega=parfourier(s12f_tao);
%     temperature_list=[temperature_list temperature];
    
    
    p1=log(-D./(omega_n.^2));
    p1=sum(p1(1*lambda+1:3*lambda));
    p2=s11_omega.*g11_omega-s12_omega.*g12_omega;
    p2=sum(p2(1*lambda+1:3*lambda));
    free_eng_list=[free_eng_list -temperature*real(2*log(2)+p1+1.5*p2)];
    temperature_list=[temperature_list temperature];
    temperature=temperature-testep
    x=0.5;
end
% 
temperature=low_temp
% 
while temperature<high_temp
    dif=Inf;
    beta=1/temperature;
    wstep = pi/points;
    w1 = zeros(1,2*points-1);
    w1(1:points)=wstep*(1:points);
    w1(points:2*points-1) = wstep*((points:2*points-1)-2*points);
    omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
    iter=0;
    tstep = beta/(2*lambda);
    while dif>err&&iter<5000
        %g_omega=parfourier(g_tao);
        s11_tao=-J.*J.*g11_tao.^2.*fliplr(g11_tao);
        s12_tao=J.*J.*g12_tao.^2.*fliplr(g12_tao);
        
        s11_omega=parfourier(s11_tao);
        s12_omega=parfourier(s12_tao);
        D=(-(exp(1i.*omega_n*tstep)-1)/tstep-mu-s11_omega).^2+(1i.*kappa-s12_omega).^2;
        
        g11_omega=(-(exp(1i.*omega_n*tstep)-1)/tstep-mu-s11_omega)./D;
        gll_tao_new=(1-x).*g11_tao+x.*parinvfourier(g11_omega);
        
        g12_omega=(-1i.*kappa+s12_omega)./D;
        gl2_tao_new=(1-x).*g12_tao+x.*parinvfourier(g12_omega);
        nd=gather(sum(abs(gll_tao_new-g11_tao).^2));

        if nd<dif
            dif=nd;
            g11_tao=real(gll_tao_new);
            g12_tao=1i*imag(gl2_tao_new);
            plot(g11_tao);
            pause(0.1);
            iter=iter+1;
        else
            x=x/2;
            continue 
        end
    end
     D=(-1i.*omega_n-mu-s11_omega).^2+(1i.*kappa-s12_omega).^2;
%     D=(-1i.*omega_n-mu-s11_omega).^2+(1i.*kappa-s12_omega).^2;
%     s11f_omega=parfourier(s11_tao);
%     g11f_omega=(-1i.*omega_n-mu-s11f_omega)./D;
%     g11f_tao=g11_tao;
%     s11f_tao=parinvfourier(s11f_omega);
%     
% 
%     g12f_omega=parfourier(s12_tao);
%     g12f_tao=g12_tao;
%     s12f_tao=s12_tao;
%     s12f_omega=parfourier(s12f_tao);
%     temperature_list=[temperature_list temperature];
    
    
    p1=log(-D./(omega_n.^2));
    p1=sum(p1(1*lambda+1:3*lambda));
    p2=s11_omega.*g11_omega-s12_omega.*g12_omega;
    p2=sum(p2(1*lambda+1:3*lambda));
    free_eng_list=[free_eng_list -temperature*real(2*log(2)+p1+1.5*p2)];
    temperature_list=[temperature_list temperature];
    temperature=temperature+testep
    x=0.5;
end
hold off
plot(temperature_list,free_eng_list)
tt=toc