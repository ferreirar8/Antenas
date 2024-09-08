clear;
clf;
close all;

% VALUES FOR THE COMPUTATION
V=1; %input('Input voltage = ');
lambda=1; %input('Wavelength = ');
a=0.005*lambda; %input('Dipole radius = ');
L=0.47*lambda; %input('Dipole length = ');
beta=(pi/2); %input('Angle of the corner refletor = '); this input doesn't work because the kernel is being calculated for N=2
S=0.5*lambda; %input('Distance between the corner and dipole = ');
k0=2*pi/lambda;
z0=120*pi;
N=pi/beta;
Ns=51;

%computation of ro
ro=zeros(1,2*N);
for i=1:2*N
    if i==1
        ro(i)=a;
    else
        ro(i)=S*sqrt(power(cos(i*(pi/2)-(pi/2))-1,2)+power(sin(i*(pi/2)-(pi/2)),2));
    end
end


delta=L/(2*Ns); %computation of delta

%computation of zm
zm=zeros(1,Ns);
for n=1:Ns
  zm(n) = (n-0.5)*delta;
end

%computation of the Z matrix
Z=zeros(Ns,Ns);
for n=1:Ns
    K=@(z) (1/(4*pi))*...                                                           % Cálculo
        ((exp(-1i*k0*sqrt(ro(1)^2+(zm(n)-z)^2))/sqrt(ro(1)^2+(zm(n)-z)^2))...       % de K
        -(exp(-1i*k0*sqrt(ro(2)^2+(zm(n)-z)^2))/sqrt(ro(2)^2+(zm(n)-z)^2))...
        +(exp(-1i*k0*sqrt(ro(3)^2+(zm(n)-z)^2))/sqrt(ro(3)^2+(zm(n)-z)^2))...
        -(exp(-1i*k0*sqrt(ro(4)^2+(zm(n)-z)^2))/sqrt(ro(4)^2+(zm(n)-z)^2)))...
        +(1/(4*pi))*...
        ((exp(-1i*k0*sqrt(ro(1)^2+(zm(n)+z)^2))/sqrt(ro(1)^2+(zm(n)+z)^2))...
        -(exp(-1i*k0*sqrt(ro(2)^2+(zm(n)+z)^2))/sqrt(ro(2)^2+(zm(n)+z)^2))...
        +(exp(-1i*k0*sqrt(ro(3)^2+(zm(n)+z)^2))/sqrt(ro(3)^2+(zm(n)+z)^2))...
        -(exp(-1i*k0*sqrt(ro(4)^2+(zm(n)+z)^2))/sqrt(ro(4)^2+(zm(n)+z)^2)));
    for m=1:Ns
        Z(n,m)=integral(K,zm(m)-(delta/2),zm(m)+(delta/2),'ArrayValued',true);  %preencher colunas com os integrais de K
        if m==Ns
            Z(n,m)=-cos(k0*zm(n));      %colocar cossenos na ultima coluna
        end
    end
end


%Definition of the sine vector
senos=zeros(1,Ns);
for m=1:Ns
    senos(m)=1i*V/(2*z0)*sin(k0*abs(zm(m)));
end
senos=senos';

IC=linsolve(Z,senos);

IC(end,:) = [];    %retirar a constante


for i=0:2
    S=0.5+i*0.5;
    disp(S);
    phi=linspace(-pi/4,pi/4,100);
    teta=pi/2;
    x=linspace(0,4);
    y=linspace(0,4);
    F_n = 2*(cos(k0.*S.*cos(phi).*sin(teta))-cos(k0.*S.*sin(phi).*sin(teta)));
    F_n = abs(F_n);
    [x_F_n, y_F_n] = pol2cart(phi, F_n);
    figure(i+1);
    plot(x_F_n, y_F_n,'LineWidth',2);
    xlim([-4 4]);
    ylim([-4 4]);
    grid on;
    hold on;
    plot(x,y,'r--');
    plot(x,-y,'r--');
    hold on;
    c = [-4 4];
    d = [0 0];
    line(c,d,'Color','black','LineWidth',2)
    c = [0 0];
    d = [-4 4];
    line(c,d,'Color','black','LineWidth',2)
    xlabel('x','fontweight','bold','fontsize',16);
    ylabel('y','fontweight','bold','fontsize',16);
    title(['Radiation pattern of the corner reflector in xoy plane with S=' num2str(S)]);
end


for i=0:2
    S=0.5+i*0.5;
    disp(S);
    phi=0;
    teta=linspace(0,pi,100);
    x=linspace(0,4);
    y=linspace(0,4);
    F_n = 2*(cos(k0.*S.*cos(phi).*sin(teta))-cos(k0.*S.*sin(phi).*sin(teta)));
    F_n = abs(F_n);
    [x_F_n, y_F_n] = pol2cart(teta, F_n);
    figure(i+4);
    plot(y_F_n, x_F_n,'LineWidth',2);
    xlim([-4 4]);
    ylim([-4 4]);
    grid on;
    hold on;
    c = [-4 4];
    d = [0 0];
    line(c,d,'Color','black','LineWidth',2)
    c = [0 0];
    d = [-4 4];
    line(c,d,'Color','black','LineWidth',2)
    xlabel('x','fontweight','bold','fontsize',16);
    ylabel('z','fontweight','bold','fontsize',16);
    title(['Radiation pattern of the corner reflector in xoz plane with S=' num2str(S)]);
    
    
    
    
end

