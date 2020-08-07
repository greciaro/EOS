function [] = PR()
% 1.- Methane, 2.- n-Butane, 3.- n-Decane
%Imput data
R = 83.1439;%cm^3-bar/K-gmol
P = 4.14;%MPa
P = P*10;%bar
T = 71.1;%C
T = T+273.15;%K
zi = [.35 .45 .2];
Ki = [6 .4 .0023];
Tc = [190.6 425.2 617.7];%K
Tr = T./Tc;
Pc = [4.54 3.8 2.12];%Mpa
Pc = Pc*10;%bar
wi = [.008 .199 .489];
M = [16.04 58.12 142.29];%g/mole

%Rachford-Rice Equation 
syms x;
eqn = sum((zi.*(1-Ki))./(Ki+(1-Ki).*x)) == 0;
solx = vpasolve(eqn,x,[0,1]);
l = double(solx);
xi = zi./(l+Ki.*(1-l));
yi = Ki.*zi./(l+Ki.*(1-l));

%PENG-ROBINSON (PR)  EOS
u = 2;
w = -1;
k12 = .0133;
k21 = k12;
k13 = .0422;
k31 = k13;
k23 = .0078;
k32 = k23;
kij = [0 k12 k13; k21 0 k23; k31 k32 0];
bi = .0778*R.*Tc./Pc;
fwi = .37464+1.54226.*wi-.26992.*(wi.^2);
ai = .45724*(R^2).*(Tc.^2)./Pc;
ai = ai.*(1+fwi.*(1-Tr.^.5)).^2;

% For Liquid
for i=1:length(ai);
aL(i) = sum(xi(i).*xi.*((ai(i).*ai).^.5).*(1-kij(i,:)));
end
aL =sum(aL);
bL = sum(xi.*bi);
AL = aL*P/((R^2)*(T^2));
BL = bL*P/(R*T);
a = 1;
b = -(1+BL-u*BL);
c = AL+w*BL^2-u*BL-u*BL^2;
d = -AL*BL-w*BL^2-w*BL^3;
EqL = [a b c d];
solnPRL = roots(EqL);
ZL = solnPRL(3);
%Volume
VL = ZL*R*T/P;
