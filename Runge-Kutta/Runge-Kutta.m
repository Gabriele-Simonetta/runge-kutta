clc
clear all
close all
%% Variabili
n1=4; %numero luci di fondo cassa 1
b1=5; % [m] larghezza luce di fondo cassa 1
a1=2.5; % [m] altezza luce di fondo cassa 1
zf1=37.5; %[m]Quota di base della luce di fondo cassa 1
L1= 10; %[m] Lunghezza dello scarico di fondo cassa 1
DZ= 0.5; % [m] Differenza quote ingresso uscita scarico
g=9.806; %[m/s²] accellerazione di gravità
chi=0; % coef di efflusso
Kstri= 50;% [m^(1/3)/s] Coefficiente di liscezza di Strickler
Ls= 150; %[m] Lunghezza dello sfioratore di superficie
Zc= 46.5; %[m s.l.m] Quota del ciglio sfiorante
%stramazzo fra cassa 1 e 2
Z12=46; %[m s.l.m] Quota del ciglio sfiorante fra cassa 1 e 2
hd=1.5; %[m] Carico di progetto sullo sfioratore
L12=10; %[m] Lunghezza sfioratore fra cassa 1 e 2
c=[[0;0.5;0.7;0.8;0.9;1],[1;0.98;0.96;0.94;0.9;0.8]]; %coeff di riduzione dello stramazz
%Per Runge_kutta
time=72; %finestra osservazione 
m=1; %minuti di discretizzazione
dt=m/60; %passo della discretizzazione
%Primi calcoli
A1=a1*b1; % [m²] Area luce di fondo
Rh1=A1/(2*(a1+b1)); %[m] Raggio idraulico della luce di fondo
Coef_Chezy= Kstri*Rh1^(1/6); %Coefficiente di Chezy c (m1/2/s)
perdi_distrib1= 2*g*L1/(Coef_Chezy^2*Rh1); %Coefficiente per le perdite di carico distribuite = 2gL/Rc2
T=load('tempi_ritorno.txt'); %matrice idrogrammi sintetici 
dat=load('Dati_casse_ZSV.txt');
Z=dat(:,1); %quota m s.l.m
V1=dat(:,3); %Volumi cassa 1
S1=dat(:,2); %superficie disponibile in cassa 1
V2=dat(:,5); %Volumi cassa 1
S2=dat(:,4); %superficie disponibile in cassa 1
Temporitorno=[T(1,:)]';
Temporitorno=Temporitorno(2:length(Temporitorno)); %facciamo ordine nel definire quali sono i tempi di ritorno
T=T(2:length(T),:);
%% Sistemiamo dove mancano i dati
TF = isnan(S1);
for i=1:length(S1)-1
    if TF(i)==1
       S1(i)=(S1(i-1)+S1(i+1))/2;
    else
       S1(i)=S1(i);
    end
    
    i=i+1;
end
if TF(length(S1))==1
    S1(length(S1))=((S1(length(S1)-1)-S1(length(S1)-2)))+S1(length(S1)-1);
end
TF = isnan(V1);
for i=1:length(V1)-1
    if TF(i)==1
       V1(i)=(V1(i-1)+V1(i+1))/2;
    else
       V1(i)=V1(i);
    end
    i=i+1;
end
if TF(length(V1))==1
    V1(length(V1))=((V1(length(V1)-1)-V1(length(V1)-2)))+V1(length(V1)-1);
end
TF = isnan(S2);
for i=1:length(S2)-1
    if TF(i)==1
       S2(i)=(S2(i-1)+S2(i+1))/2;
    else
       S2(i)=S2(i);
    end
    
    i=i+1;
end
if TF(length(S2))==1
    S2(length(S2))=((S2(length(S2)-1)-S2(length(S2)-2)))+S2(length(S2)-1);
end
TF = isnan(V2);
for i=1:length(V2)-1
    if TF(i)==1
       V2(i)=(V2(i-1)+V2(i+1))/2;
    else
       V2(i)=V2(i);
    end
    
    i=i+1;
end
if TF(length(V2))==1
    V2(length(V2))=((V2(length(V2)-1)-V2(length(V2)-2)))+V2(length(V2)-1);
end

%% Tipo di funzionamento
adisc1=a1*1.5+zf1; %discriminante
for i=1:length(Z)
    if Z(i)<adisc1
        Qout(i)=(2/((3+chi)^(3/2)))*4*b1*(Z(i)-zf1)*(2*g*(Z(i)-zf1))^(1/2);
    elseif Z(i)<Zc
        Qout(i)=(1/((1+chi+perdi_distrib1)^(1/2)))*4*b1*a1*(2*g*((Z(i)-zf1)+DZ-a1))^(1/2);
    else
        csi(i)=(Z(i)-Zc)/hd;
        mu(i)=(2/(3*(3^(1/2))))*(1+(4*csi(i))/(9+5*csi(i)));
        Qout(i)=((1/((1+chi+perdi_distrib1)^(1/2)))*4*b1*a1*(2*g*((Z(i)-zf1)+DZ-a1))^(1/2))+ mu(i)*Ls*(Z(i)-Zc)*(2*g*(Z(i)-Zc))^(1/2);
    end
end
Qout=Qout';
ZQSV=[Z,Qout,S1,V1,S2,V2];
fileID = fopen('ZQSV.txt','w');
fprintf(fileID,'%12.4f%12.4f %12.4f%12.4f%12.4f %12.4f\n',ZQSV');
fclose(fileID);  
%% Metodo Runge - Kutta
%condizioni iniziali
quota1(1) = Z(1);   
quota2(1) = Z(5);
% for j=1:7 % ciclo per ciascun tempo di ritorno
% tempo=j;
% %Quale tempo di ritorno?
% if tempo==1
%     QT=T(:,2);
% elseif tempo==2
%     QT=T(:,3);
%     elseif tempo==3
%     QT=T(:,4);
%     elseif tempo==4
%     QT=T(:,5);
%     elseif tempo==5
%     QT=T(:,6);
%     elseif tempo==6
%     QT=T(:,7);
%     elseif tempo==7
%     QT=T(:,8);
% end
Qe=[T(:,2)];    %Qe=[QT];
Q12(1)=0;
Sup1(1)=S1(2);
Sup2(1)=S2(5);
%% Discretizzazione Qe per tempi (t+Delta t)
Qdt=Qe*dt;
i=1;
for t=dt:dt:time %!! VA DA 0 A time o da dt a time?
  if i ==1
      Qedt(i)= Qdt(1);
    elseif i>1
       Qedt(i)= Qdt(ceil(dt))+Qedt(i-1);
      end
 i=i+1;
end
%Discretizzazione Qe per tempi (t+Delta t)
Qdtmezz=Qe*dt/2;
i=1;
for t=dt/2:dt/2:time %!! VA DA 0 A time o da dt a time?
  if i ==1
      Qedtmezz(i)= Qdtmezz(1);
    elseif i>1
       Qedtmezz(i)= Qdtmezz(ceil(dt/2))+Qedtmezz(i-1);
      end
 i=i+1;
end
%% Runge_Kutta
tt=dt:dt:time; %appoggio per interpolazioni
ttmez=dt/2:dt/2:time;%appoggio per interpolazioni
i=1;% i usato come contatore step
for t=dt:dt:time
    
%% PER IL CODICE COMPLETO CONTATTAMI
