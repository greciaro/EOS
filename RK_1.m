function [RK] = RK()
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


%REDLICH-KWONG (RK)EOS
u = 1;
w = 0;
ai = .42748*(R^2).*(Tc.^(5/2))./(Pc.*(T^(1/2)));
bi = .08664*R.*Tc./Pc;

% For Liquid
for i=1:length(ai);
aL(i) = sum(xi(i).*xi.*((ai(i).*ai).^(1/2)));
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
solnRKL = roots(EqL);
ZL = solnRKL(3);
%Volume
VL = ZL*R*T/P;

%Fugacity
deltaLi = 2.*(ai./aL).^.5;
uu = (u^2-4*w)^.5;
bi_bL = (Tc./Pc)./sum(xi.*Tc./Pc);
coefugLi = bi_bL.*(ZL-1)-log(ZL-BL);
coefugLi1 = AL/(BL*(u^2-4*w)^.5).*(bi_bL-deltaLi);
coefugLi2 = 2*ZL+BL*(u+uu);
coefugLi3 = 2*ZL+BL*(u-uu);
coefugLi = coefugLi + coefugLi1.*log(coefugLi2/coefugLi3);
coefugLi = exp(coefugLi);
fugLi = P.*coefugLi.*xi;

