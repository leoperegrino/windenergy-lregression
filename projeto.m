
%%%%dados referentes ao ecmwf%%%

ecmwf=ncgeodataset('dados_ecmwf.nc');

u=ecmwf.data{'10_metre_U_wind_component_surface'};  %carregar velocidade u
u(:,:,1)=[];                                        %retirar dados da longitude mais afastada
u=double(u);
v=ecmwf.data{'10_metre_V_wind_component_surface'};  %carregar velocidade v
v(:,:,1)=[];                                        %retirar dados da longitude mais afastada
v=double(v);

%variaveis em formato timetable

tempo=datetime(datevec(ecmwf.time('time')),'TimeZone','+00:00');
gcm=timetable(tempo,u,v);
gcm.tempo.TimeZone='-03:00';

clear u v tempo ecmwf







%%%dados referentes ao anemômetro%%%

v=load('T.mat');           %carregar velocidades
v=v.T;                     %atribuir velocidades a uma matriz

%tempo em formato datetime

tempo=load('DATA_torre.mat');
tempo=datetime(datevec(tempo.DATA_torre,'yyyymmddHHMM'),'TimeZone','-03:00');

anem=timetable(tempo,v);
%anen.v=anem.v(~isnan(anem.v));
clear tempo v 



% anemn=timetable(anem.tempo(deltat_anem),anem.v(deltat_anem));
% anemm=timetable(gcm.tempo(deltat_ecmwf),y_reg);


%%%regressão linear múltipla%%%

m=hypot(gcm.u,gcm.v);                                %tirar módulo de (u,v)
deltat_anem=(28:36:length(anem.v)).';                %criar contador de 6 em 6h
deltat_ecmwf=(11693:11692+length(deltat_anem)).';    %contador T0=T0_anem


M=m(deltat_ecmwf,:,:);
Y=anem.v(deltat_anem);     %definir vetor de variáveis e coef linear
X=[ones(length(deltat_anem),1) M(:,1,1) M(:,1,2) M(:,2,1) M(:,2,2)];


[coef,~,~,~,~]=regress(Y,X);      	  %]regressão do 
reg_anem=fitlm(X(:,2:5),Y);           %]anem

disp(reg_anem.Rsquared);
clear deltat_anem deltat_ecmwf
keyboard %pausa na operação

y_reg=coef(1)+coef(2)*M(:,1,1)+coef(3)*M(:,1,2)+coef(4)*M(:,2,1)+coef(5)*M(:,2,2);

scatter(y_reg,Y)
ylabel({'Velocidade do anem','(m/s)'});
xlabel({'Velocidade de regressão do anem','(m/s)'});
lsline;
legend({'$V_{anem}(\hat{Y}_{anem})$','Ajuste Linear'},'Interpreter','latex','Location','southeast','FontSize',10);

v_final=[y_reg;Y];                               %conjunto de velocidades anem+ecmwf

wb=fitdist(v_final,'Weibull');


keyboard %pausa na operação
clear deltat_anem deltat_ecmwf X coef m y_reg M Y gcm anem







%%%cálculo do raio e Udesign%%%


u=0:0.5:25;                                                                                                              %]calcular EP 
dens_P_mec=[u; (1/2)*1.22*u(1:51).^3].';                                                                                 %|normalizada por R^2
dens_PP_mec=[u(1:50)+0.25; ((cdf(wb,u(2:51))-cdf(wb,u(1:50))).*((dens_P_mec(1:50,2)+dens_P_mec(2:51,2)).'*(1/2)))].';    %]de cada bin
                                                                                         


[~,I]=max(dens_PP_mec(:,2));                 %]
U_design=dens_PP_mec(I,1);                   %]descobrir max potência normalizada e sua velocidade (Udesign) 
                            
clear I dens_PP_mec



syms R;                                             %]
eq=(1/2)*1.22*pi*R^2*U_design^3==400000*0.35/0.44;  %|calcular raio do rotor
R=double(solve(eq,R));                              %|
R=R(R>0);                                           %]

clear eq
                        






%%%cálculo das propriedades da pá%%%


Cl=0.7725;                                                                          %]
r=(R/40:R/40:R).';                                                                   %|
W=(35*2*pi/60);                                                                     %|calcular propriedades

phi=[r (2/3)*atan((1./r)*(U_design)/W)];                                            %|da pá com Cd=0 e 
c=[r ((1/(3*Cl))*8*pi*r.*(1-cos(phi(:,2))))];                                       %|alfa(Cl/Cd(max))
a=[r 1./(1+(((2*pi*r).*(4*sin(phi(:,2)).^2)./(3*Cl*c(:,2).*cos(phi(:,2))))))];      %|
a_=[r ((1-3*a(:,2))./(4*a(:,2)-1))];                                                %]
beta=[r phi(:,2)-0.0698];

keyboard %pausa na operação










%%%curva de potência%%%


u_rel=(u(1:40).'.*(1-a(:,2))).^2+(W*r.*(1+a_(:,2))).^2;

P=ones(40,1);
for i=1:40
    P(i)=sum(Cl*sin(phi(:,2)).*c(:,2).*(r-R/80).*u_rel(i)*W*3*1.22*0.5*R/40);
end
P=[u(1:40).' P/1000];



keyboard %pausa na operação
clear Cl i W u_rel











%%%gráficos%%%


hold on
yyaxis left
plot(u(1:50)+0.25,dens_P_mec(1:50,2)*(pi*R^2)/1000,'.:');
plot(u(1:40),P(:,2),'-.');
scatter(U_design, 0.5*1.22*pi*R^2*U_design^3*0.001*0.44,'ok');
ylim([0 420]);
ylabel({'Potência','(kW)'});

yyaxis right
plot(0:.001:25,pdf(wb,0:.001:25));
ylabel({'Densidade de probabilidade','(s/m)'});
xlabel({'Velocidade','(m/s)'});
xlim([0 u(end)]);

legend({'P_{disp}','P_{mec}','P_{e,des}','p'},'Location','southeast');
hold off

keyboard %pausa na operação

plot(r,rad2deg(phi(:,2)),'-o');
ylabel({'\phi','(\circ)'});
xlabel({'Raio','(m)'});
legend({'\phi(r)'},'Location','southwest');
xlim([0 R]);

keyboard %pausa na operação

plot(r,rad2deg(beta(:,2)),'-o');
ylabel({'\beta','(\circ)'});
xlabel({'Raio','(m)'});
xlim([0 R]);
legend({'\beta(r)'},'Location','southwest');

keyboard %pausa na operação

c0=c(end,2);
hold on
plot(r,c(:,2),'-o');
fplot(@(x) c0,'-.k')
ylabel({'Corda da pá','(m)'});
xlabel({'Raio','(m)'});
xlim([0 R]);
ylim([0 4.5]);
legend({'c(r)'},'Location','east');
hold off

keyboard %pausa na operação

% plot(u(1:50)+0.25,EP_R_e(:,2)/1000,'-o');
% ylabel({'Densidade energética','(kW\timess/m)'});
% xlabel({'Velocidade','(m/s)'});
% xlim([0 u(end)]);
% legend({'p_{e}(v)'},'Location','east');
% 
% keyboard %pausa na operação

% plot(0:.001:25,pdf(wb,0:.001:25));
% ylabel({'Densidade de probabilidade','(s/m)'});
% xlabel({'Velocidade','(m/s)'});
% legend({'p(v)'},'Location','east');
% xlim([0 u(end)]);

keyboard %pausa na operação
clf

hold on
plot(r,a(:,2),'-o');
plot(r,a_(:,2),'-o');
ylim([0 0.5]);
xlim([0 R]);
ylabel({'Coeficientes de indução axial e radial','(o)'});
xlabel({'Raio','(m)'});
legend({'axial','radial'},'Location','east');
hold off

keyboard
clall


