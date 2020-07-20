clear;
tic
%p=parpool;

global beta;
global lambda;

high_temp=0.04;
low_temp=0.04;
testep=0.005;
temperature=high_temp;
err=10^(-25);

lambda=2^16; 
beta=1/temperature;
j=1;
q=4;
x=0.5;
max=2*lambda;
points=max;
mu=0.75;
omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
tstep = beta/(2*lambda);
g_tao = -0.5*sign(points+0.5-(1:2*points));
glr_tao = -0.5*sign(points+0.5-(1:2*points));
g_omega = -conj(0.5*tstep*fft(g_tao));
glr_omega = -conj(0.5*tstep*fft(g_tao));
s_tao = zeros(1,2*points);
s_omega = zeros(1,2*points);
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
        s_tao=j.*j.*(2.*g_tao).^(q-1)./q;
        s_omega=fftfourier(s_tao);
        g_omega(2:2*points)=((-(exp(-1i*w1)-1)/tstep+s_omega(2:2*points))./((-(exp(-1i*w1)-1)/tstep+s_omega(2:2*points)).^2-mu ^2)).*mod(1:4*lambda-1,2);
        g_tao_new=(1-x).*g_tao+x.*fftinvfourier(g_omega);
        nd=gather(sum(abs(g_tao_new-g_tao).^2));

        if nd<dif
            dif=nd;
            g_tao=real(g_tao_new);
            %plot(g_tao);
            %pause(0.01);
            iter=iter+1;
        else
            x=x/2;
            continue 
        end
    end
    glr_omega(2:2*points)=((-1i*mu)./(((-(exp(-1i*w1)-1)/tstep+s_omega(2:2*points)).^2-mu ^2))).*mod(1:4*lambda-1,2);
    glr_tao=fftinvfourier(glr_omega);
    sf_omega=parfourier(s_tao);
    gf_omega=(-(1i.*omega_n+sf_omega)./((1i.*omega_n+sf_omega).^2-mu^2));
    temperature_list=[temperature_list temperature];
    p1=log(((1+sf_omega./(1i.*omega_n)).^2+mu.^2./(omega_n).^2));
    p1=gather(sum(p1(lambda+1:3*lambda+1)));
    p2=real(gather(sum(sf_omega(lambda+1:3*lambda+1).*gf_omega(lambda+1:3*lambda+1))));
    p3=gather(sum(g_tao(2*lambda+1:4*lambda).^q))
    free_eng_list=[free_eng_list real(-log(2)./(2.*beta)-1/(4*beta)*p1-1/(2.*beta)*p2-(beta*j^2*2^(q-1)*p3)/(4*lambda*q^2))];
    temperature=temperature-testep
    x=0.5;
end
% 
% temperature=low_temp
% 
% while temperature<high_temp
%     dif=Inf;
%     beta=1/temperature;
%     
%     wstep = pi/points;
%     w1 = zeros(1,2*points-1);
%     w1(1:points)=wstep*(1:points);
%     w1(points:2*points-1) = wstep*((points:2*points-1)-2*points);
%     omega_n= (2.*pi.*((1:2*max) - max - 0.5))./beta;
%     iter=0;
%     tstep = beta/(2*lambda);
%     while dif>err&&iter<5000
%         %g_omega=parfourier(g_tao);
%         s_tao=j.*j.*(2.*g_tao).^(q-1)./q;
%         s_omega=fftfourier(s_tao);
%         g_omega(2:2*points)=((-(exp(-1i*w1)-1)/tstep+s_omega(2:2*points))./((-(exp(-1i*w1)-1)/tstep+s_omega(2:2*points)).^2-mu^2)).*mod(1:4*lambda-1,2);
%         g_tao_new=(1-x).*g_tao+x.*fftinvfourier(g_omega);
%         nd=gather(sum(abs(g_tao_new-g_tao).^2));
% 
%         if nd<dif
%             dif=nd;
%             g_tao=real(g_tao_new);
%             plot(g_tao);
%             pause(0.01);
%             iter=iter+1;
%         else
%             x=x/2;
%             continue
%         end
%     end
%     sf_omega=parfourier(s_tao);
%     gf_omega=(-(1i.*omega_n+sf_omega)./((1i.*omega_n+sf_omega).^2-mu^2));
%     temperature_list=[temperature_list temperature];
%     p1=log(((1+sf_omega./(1i.*omega_n)).^2+mu.^2./(omega_n).^2));
%     p1=gather(sum(p1(lambda+1:3*lambda+1)));
%     p2=real(gather(sum(sf_omega(lambda+1:3*lambda+1).*gf_omega(lambda+1:3*lambda+1))));
%     p3=gather(sum(g_tao(2*lambda+1:4*lambda).^q))
%     free_eng_list=[free_eng_list real(-log(2)./(2.*beta)-1/(4*beta)*p1-1/(2.*beta)*p2-(beta*j^2*2^(q-1)*p3)/(4*lambda*q^2))];
%     temperature=temperature+testep
%     x=0.5;
% end
hold off
plot(g_tao);
hold on 
plot(imag(glr_tao));
legend('g_tao','imag(glr_tao)')
tt=toc