function [Q12]=stram12(Z1,Z2,Z12,c,hd,L12,g)
csi2=(Z1-Z12)/hd;
if csi2<0
    Q12=0;
else
    yi=(Z2-Z12)/(Z1-Z12);
    ci=interp1(c(:,1),c(:,2),yi,'linear'); %coefficiente di rigurgito
    mu12=(2/(3*(3^(1/2))))*(1+(4*csi2)/(9+5*csi2));
    Q12= ci*mu12*L12*(Z1-Z12)*(2*g*(Z1-Z12))^(1/2); %portata fra cassa 1 e 2
end
end