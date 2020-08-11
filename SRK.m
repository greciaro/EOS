function [] = SRK()
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

%SOAVE-REDLICH-KOWNG (SRK) EOS
u = 1;
w = 0;
k12 = .0056;
k13 = .0411;
k23 = .0067;
k21 = k12;
k31 = k13;
k32 = k23;
kij = [0 k12 k13; k21 0 k23; k31 k32 0];
fwi = .48+1.574.*wi-.176.*(wi.^2);
ai = .42748*(R^2).*(Tc.^2)./Pc;
ai = ai.*((1+fwi.*(1-(Tr.^(1/2)))).^2);
bi = .08664*R.*Tc./Pc;

% For Liquid
for i=1:length(ai);
aL(i) = sum(xi(i).*xi.*((ai(i).*ai).^(1/2)).*(1-kij(i,:)));
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
solnSRKL = roots(EqL);
ZL = solnSRKL(3);
%Volume
VL = ZL*R*T/P;

%Fugacity
for i=1:length(ai);
deltaLi1(i) = sum(xi.*(ai.^.5).*(1-kij(i,:)));
deltaLi(i) = (2.*(ai(i).^.5)./aL)*deltaLi1(i);
end
uu = (u^2-4*w)^.5;
bi_bL = (Tc./Pc)./sum(xi.*Tc./Pc);
coefugLi = bi_bL.*(ZL-1)-log(ZL-BL);
coefugLi1 = AL/(BL*(u^2-4*w)^.5).*(bi_bL-deltaLi);
coefugLi2 = 2*ZL+BL*(u+uu);
coefugLi3 = 2*ZL+BL*(u-uu);
coefugLi = coefugLi + coefugLi1.*log(coefugLi2/coefugLi3);
coefugLi = exp(coefugLi);
fugLi = P.*coefugLi.*xi;

% For Vapor
for i=1:length(ai);
aV(i) = sum(yi(i).*yi.*((ai(i).*ai).^(1/2)).*(1-kij(i,:)));
end
aV =sum(aV);
bV = sum(yi.*bi);
AV = aV*P/((R^2)*(T^2));
BV = bV*P/(R*T);
a = 1;
b = -(1+BV-u*BV);
c = AV+w*BV^2-u*BV-u*BV^2;
d = -AV*BV-w*BV^2-w*BV^3;
EqV = [a b c d];
solnSRKV = roots(EqV);
ZV = solnSRKV(1);
%Volume
VV = ZV*R*T/P;

%Fugacity
for i=1:length(ai);
deltaVi1(i) = sum(yi.*(ai.^.5).*(1-kij(i,:)));
deltaVi(i) = (2.*(ai(i).^.5)./aV)*deltaVi1(i);
end
bi_bV = (Tc./Pc)./sum(yi.*Tc./Pc);
coefugVi = bi_bV.*(ZV-1)-log(ZV-BV);
coefugVi1 = AV/(BV*(u^2-4*w)^.5).*(bi_bV-deltaVi);
coefugVi2 = 2*ZV+BV*(u+uu);
coefugVi3 = 2*ZV+BV*(u-uu);
coefugVi = coefugVi + coefugVi1.*log(coefugVi2/coefugVi3);
coefugVi = exp(coefugVi);
fugVi = P.*coefugVi.*yi;

%Results
fprintf('\n');
fprintf('**********************************');
fprintf('\n');
fprintf('      SOAVE-REDLICH-KOWNG EOS');
fprintf('\n');
fprintf('\n');
fprintf(' LIQUID ');
fprintf('\n');
fprintf('Volume = %d',VL);
fprintf('\n');
fprintf('f_methane = %d',fugLi(1));
fprintf('\n');
fprintf('f_n-butaneane = %d',fugLi(2));
fprintf('\n');
fprintf('f_n-decane = %d',fugLi(3));
fprintf('\n');
fprintf('\n');
fprintf(' VAPOR ');
fprintf('\n');
fprintf('Volume = %d',VV);
fprintf('\n');
fprintf('f_methane = %d',fugVi(1));
fprintf('\n');
fprintf('f_n-butaneane = %d',fugVi(2));
fprintf('\n');
fprintf('f_n-decane = %d',fugVi(3));
fprintf('\n');
fprintf('\n');
fprintf('**********************************');
fprintf('\n');
fprintf('\n');

end

