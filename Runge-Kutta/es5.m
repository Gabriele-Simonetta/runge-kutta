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
    
 if i==1
    Q1(1)=Qedt(1);
    DZ11=dt*(Qedt(2)-Q1(1)-Q12(1))/(Sup1(1)); %!! Qedt(2)  cioè allo step successivo altrimenti  rimangono nulli
    DZ21=dt*(Q12(1))/(Sup2(1));
 else 
    Sup1=interp1(Z,S1,Z1,'linear');%S1 alla quota Z1
    Qout1=interp1(Z,Qout,Z1,'linear');%Q1 alla quota Z1
    Qing=interp1(ttmez,Qedtmezz,(t+dt/2),'linear');
    [Q12]=stram12(Z1,Z2,Z12,c,hd,L12,g);
    DZ11=dt*(Qing-Qout1-Q12)/(Sup1);
    DZ21=dt*(Q12-0)/(Sup2);
  end

 %secondo passo
 Z1=quota1+DZ11/2;
 Sup1=interp1(Z,S1,Z1,'linear');%S1 alla quota Z1
 Qout1=interp1(Z,Qout,Z1,'linear');%Q1 alla quota Z1
 Qing=interp1(ttmez,Qedtmezz,(t+dt/2),'linear'); %Qe alla quota Z1
 
Z2=quota2+DZ21/2; %!!Possibile errore fra quota pelo libero e quota invaso
%Calcolo Q12 
[Q12]=stram12(Z1,Z2,Z12,c,hd,L12,g);
DZ12=dt*(Qing-Qout1-Q12/Sup1);
DZ22=dt*(Q12-0/Sup2);

%%terzo passo

%Soluzione
%quota1(i)=quota1(i-1)+(DZ11+2*DZ12+2*DZ13+DZ14)/6;
%quota2(i)=quota2(i-1)+(DZ21+2*DZ22+2*DZ23+DZ24)/6;
i=i+1;
end
    
%non  cancellare più in basso sono già impostate le immagini che si
%dovranno salvare


% QT=T(:,2)
% for i=1:length(T)
% DZ1(i)=dt*(QT-Qout(i)
% end
% end
% 


% W=W'; %volume invasabil
% % t=T(:,1);
% % dt=3600/2;
% % O=W*10^3+Qu*dt;
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
% for i=1:length(T)
%     if i==1
%         Mat(1,4,j)=QT(1,1);
%         zi = interp1(Qu,Z,QT(1,1),'linear');
%         Wi = interp1(Qu,W,QT(1,1),'linear');
%         Mat(1,2,j)=zi;
%         Mat(1,3,j)=Wi;
%         Mat(1,1,j)=0;
%     else
%         Mat(i,1,j)=dt*(QT(i,1)+QT(i-1,1))-dt*Mat(i-1,4,j)+Mat(i-1,3,j)*10^3;
%         Mat(i,2,j) = interp1(O,Z,Mat(i,1,j),'linear');
%         Mat(i,3,j) = interp1(O,W,Mat(i,1,j),'linear');
%         Mat(i,4,j) = interp1(O,Qu,Mat(i,1,j),'linear');    
%     end
%     i=i+1;
% end
% QeMAX(j)=max (QT);
% QuMAX(j)=max ( Mat(:,4,j));
% ZMAX(j)=max ( Mat(:,2,j));
% WMAX(j)=max ( Mat(:,3,j));
% effi(j)=1-(QuMAX(j)/QeMAX(j));
% ridotta(j)=-log(log(Temporitorno(j)/(Temporitorno(j)-1)));
% j=j+1;
% end
% %% salvataggio dati
% R=[Temporitorno',QeMAX',QuMAX',effi',ZMAX',WMAX'];
% fileID = fopen('Risultati.txt','w');
% fprintf(fileID,'%12.4f %12.4f %12.4f %12.4f %12.4f %12.4f\n',R');
% fclose(fileID);  
% %% Salvataggio figure
% %Curve aree-volumi
% fg=1;
% 
% figure(fg);
% line(W(1:31)*10^-3,Z(1:31),'Color','r')
% ax1 = gca; % current axes
% ax1.XColor = 'r';
% ax1.YColor = 'k';
% ax1_pos = ax1.Position; % position of first axes
% xlabel ('Volume invasabile [10^6 m³]');
% ylabel ('Quota idrica [m s.l.m]');
% grid on
% ax2 = axes('Position',ax1_pos,'XAxisLocation','top','YAxisLocation','right','Color','none');
% set(gca,'ytick',[])%delete second y-axis
% xlabel ('Aree [ha]');
% line(S1(1:31),Z(1:31),'Parent',ax2,'Color','k');
% ax2.XTick = 0:10:150;
% grid on
% saveas(gcf,'Curve_aree-volumi_fig.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'Curve_aree-volumi_fig.png'); %salvataggio grafico formato .png
% 
% %Interpolazione a mano
% fg=fg+1;
% figure(fg);
% plot(O(1:31),Z(1:31));
% grid on
% grid minor
% xlabel('O [m³]');
% ylabel('Quota idrica [m s.l.m]');
% saveas(gcf,'hinterp_fig.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'hinterp_fig.png'); %salvataggio grafico formato .png
% 
% 
% fg=fg+1;
% figure(fg);
% plot(O(1:31),Qu(1:31));
% grid on
% grid minor
% xlabel('O [m³]');
% ylabel('Portata uscente [m³/s]');
% saveas(gcf,'Qinterp_fig.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'Qinterp_fig.png'); %salvataggio grafico formato .png
% 
% fg=fg+1;
% figure(fg);
% plot(O(1:31),W(1:31));
% grid on
% grid minor
% xlabel('O [m³]');
% ylabel('Volume invasabile [m³]');
% saveas(gcf,'Winterp_fig.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'Winterp_fig.png'); %salvataggio grafico formato .png
% 
% %Volumi-Quota onda di 50 anni
% fg=fg+1;
% figure(fg);
% hold on
% yyaxis left
% plot(T(:,1),Mat(:,3,3));
% title('T=50 anni');
% xlabel('Tempo [h]')
% ylabel('Volumi [10³ m³]');
% yyaxis right
% plot(T(:,1),Mat(:,2,3));
% ylabel('Quota idrica [m s.l.m]');
% legend('Volumi invasati','Quota pelo libero');
% grid on
% grid minor
% hold off
% saveas(gcf,'Volumi+Quota_50anni_fig.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'Volumi+Quota_50anni_fig.png'); %salvataggio grafico formato .png
% 
% %200 anni
% fg=fg+1;
% figure(fg);
% plot(T(:,1),T(:,7),'k',T(:,1),Mat(:,4,6),'--k');
% grid on
% grid minor
% xlabel('Tempo[h]');
% ylabel('Portata  [m³/s]');
% ylim([0 (max(QeMAX)+100)]);
% xlim([0 (i+1)]);
% legend('Portata entrante T=200','Portata uscente T=200');
% saveas(gcf,'200anni.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'200anni.png'); %salvataggio grafico formato .png
% 
% fg=fg+1;
% figure(fg);
% hold on
% yyaxis left
% plot(T(:,1),Mat(:,3,6));
% title('T=200 anni');
% xlabel('Tempo [h]')
% ylabel('Volumi [10³ m³]');
% yyaxis right
% plot(T(:,1),Mat(:,2,6));
% ylabel('Quota idrica [m s.l.m]');
% legend('Volumi invasati','Quota pelo libero');
% grid on
% grid minor
% hold off
% saveas(gcf,'Volumi+Quota_200anni_fig.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'Volumi+Quota_200anni_fig.png'); %salvataggio grafico formato .png
% 
% %Idrogrammi in ingresso
% fg=fg+1;
% figure(fg);
% plot(T(:,1),T(:,2),'y',T(:,1),T(:,3),'m',T(:,1),T(:,4),'b',T(:,1),T(:,5),'c',T(:,1),T(:,6),'r',T(:,1),T(:,7),'k',T(:,1),T(:,8));
% grid on
% grid minor
% xlabel('Tempo[h]');
% ylabel('Portata entrante [m³/s]');
% ylim([0 (max(QeMAX)+100)]);
% xlim([0 (i+1)]);
% legend((sprintf("T=5")),(sprintf("T=10")),(sprintf("T=20")),(sprintf("T=50")),(sprintf("T=100")),(sprintf("T=200")),(sprintf("T=500")));
% saveas(gcf,'Idrogrammisintetici_fig.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'Idrogrammisintetici_fig.png'); %salvataggio grafico formato .png
% 
% %idrogrammi in uscita
% fg=fg+1;
% figure(fg);
% plot(T(:,1),Mat(:,4,1),'--y',T(:,1),Mat(:,4,2),'--m',T(:,1),Mat(:,4,3),'--b',T(:,1),Mat(:,4,4),'--c',T(:,1),Mat(:,4,5),'--r',T(:,1),Mat(:,4,6),'--k',T(:,1),Mat(:,4,7),'--');
% grid on
% grid minor
% xlabel('Tempo[h]');
% ylabel('Portata uscente [m³/s]');
% ylim([0 (max(QeMAX)+100)]);
% xlim([0 (i+1)]);
% legend((sprintf("T=5")),(sprintf("T=10")),(sprintf("T=20")),(sprintf("T=50")),(sprintf("T=100")),(sprintf("T=200")),(sprintf("T=500")));
% saveas(gcf,'Ondeinuscita_fig.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'Ondeinuscita_fig.png'); %salvataggio grafico formato .png
% 
% %Quota/Portata
% fg=fg+1;
% figure(fg);
% plot(Mat(:,4,1),Mat(:,2,1),'--y',Mat(:,4,2),Mat(:,2,2),'--m',Mat(:,4,3),Mat(:,2,3),'--b',Mat(:,4,4),Mat(:,2,4),'--c',Mat(:,4,5),Mat(:,2,5),'--r',Mat(:,4,6),Mat(:,2,6),'--k',Mat(:,4,7),Mat(:,2,7),'--');
% grid on
% grid minor
% xlabel('Portate uscenti [m³/s]');
% ylabel('Quota idrica [m s.l.m]');
% legend((sprintf("T=5")),(sprintf("T=10")),(sprintf("T=20")),(sprintf("T=50")),(sprintf("T=100")),(sprintf("T=200")),(sprintf("T=500")));
% saveas(gcf,'Quota_Portata_fig.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'Quota_Portata_fig.png'); %salvataggio grafico formato .png
% 
% %Quota massima e Volume massimo invasato
% fg=fg+1;
% figure(fg);
% hold on
% yyaxis left
% plot(Temporitorno,WMAX);
% %title('Plots with Different y-Scales')
% xlabel('Tempo di ritorno')
% ylabel('Volume [10³ m³]');
% yyaxis right
% plot(Temporitorno,ZMAX);
% ylabel('Quota idrica [m s.l.m]');
% legend('Volumi invasati','Quota pelo libero');
% grid on
% grid minor
% hold off
% saveas(gcf,'Volumi+Quota_fig.fig'); %salvataggio grafico formato .fig
% saveas(gcf,'Volumi+Quota_fig.png'); %salvataggio grafico formato .png
% 
% %efficienza
% fg=fg+1;
% figure(fg)
% ymax=0.6;
% [x,y,ridotta] = carta_di_gumbel(Temporitorno,ymax);
% hold on
% plot(ridotta,effi);
% ylabel('Efficienza');
% grid 'on'
% grid minor
% xlabel ('Variabile ridotta di Gumbel');
% for i=1:length(Temporitorno)
% line(x(i,:),y,'Color','k');
% text(x(i,1),ymax,texlabel(strcat('T=',num2str(Temporitorno(i)))));
% i=i+1;
% end
% hold off
% saveas(gcf,'Efficienza_fig.fig');%salvataggio grafico formato .fig
% saveas(gcf,'Efficienza_fig.png');%salvataggio grafico formato .png
% 
% %confronto_Q
% fg=fg+1;
% figure(fg)
% ymax=1400;
% [x,y,ridotta] = carta_di_gumbel(Temporitorno,ymax);
% hold on
% plot(QeMAX,ridotta, 'y',QuMAX,ridotta,'r');
% xlabel('Portata [m³/s]');
% grid 'on'
% grid minor
% ylabel ('Variabile ridotta di Gumbel');
% for i=1:length(Temporitorno)
% line(y,x(i,:),'Color','k');
% text(ymax,x(i,1),texlabel(strcat('T=',num2str(Temporitorno(i)))));
% i=i+1;
% end
% legend('Portate in ingresso','Portate in uscita');
% hold off
% saveas(gcf,'confronto_Q_fig.fig');%salvataggio grafico formato .fig
% saveas(gcf,'confronto_Q_fig.png');%salvataggio grafico formato .png