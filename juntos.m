clear;
close all;
clf;

% VALUES FOR THE COMPUTATION
V=1; %input('Input voltage = ');
lambda=1; %input('Wavelength = ');
a=0.005*lambda; %input('Dipole radius = ');
L=0.47*lambda; %input('Dipole length = ');
beta=90; %input('Angle of the corner refletor = '); this input doesn't work because the kernel is being calculated for N=2
S=100*lambda; %input('Distance between the corner and dipole = ');
k0=2*pi/lambda;
z0=120*pi;
N=180/beta;

%computation of ro
ro=zeros(1,2*N);
for i=1:2*N
    if i==1
        ro(i)=a;
    else
        ro(i)=S*sqrt(power(cosd(i*beta-90)-1,2)+power(sind(i*beta-90),2));
    end
end

for ns=1:3
    if ns==1
        Ns=21;
    elseif ns==2
        Ns=21;
    else
        Ns=21;
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
    S=zeros(1,Ns);
    for m=1:Ns
        S(m)=-1i*V/(2*z0)*sin(k0*abs(zm(m)));
    end
    S=S';



    IC=linsolve(Z,S);

    C1 = IC(Ns,1);    %gravar a constante C1 do vetor
    IC(end,:) = [];    %retirar a constante
    
    Zin=V./IC;  %Cálculo da impedância de entrada
    
    %Fazer plot à corrente
    z=0:1:Ns;
    y=abs(IC)/(abs(IC(1)));
    y=[y;zeros(2,1)];

    figure (1)
    stairs(z/(Ns),y,'LineWidth',2)
    hold on
    title('Amplitude da corrente')
    grid
    xlabel('2z/L');
    ylabel('Amplitude I(z)');
    legend('Ns=21','Ns=51','Ns=71')
    zx=0:0.001:L/2;
    Im2t = IC(1,1);
    Im=abs(Im2t)/sin((k0*L)/2);
    Iaprox=sin(k0*((L/2)-zx));
    plot(zx/(L/2),Iaprox);
    hold on;
end

lambda=1;
S=1*lambda;
theta=0:0.05:pi;
phi=0:0.05:pi;
S=100;
DIAG=zeros(21,length(theta));
for num=1:length(theta)
    c1=exp(1i*k0*zm*cos(theta(num))); 
    c2=exp(1i*k0*(-zm)*cos(theta(num))); 
    c3=c1+c2;
    %integral
    for it=1:Ns
        soma=((-IC)/(-IC(1,1)))*c3*delta; 
    end
    %Fi0=abs(2*(cos(k0*S*cos(phi(num))*sin(theta(num)))-cos(k0*S*sin(phi(num))*sin((theta(num))))));
    he=(sin(theta(num)*soma));
    DIAG(1,num)=he(20,21);
end
disp(size(he));
disp(num);

figure (2)
polarplot(theta,abs(DIAG));
hold on;
[x,y]=pol2cart(theta,abs(DIAG));
figure (3)
plot(x,y);
hold on;
    %correntes= (IC(num)/IC(1,1))*exp(i*k0*z2*cos(theta(num)));
    
    
%polarplot(theta,he(num));