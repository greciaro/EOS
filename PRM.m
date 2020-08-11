function [] = PRM()
% 1.- Methane, 2.- n-Butane, 3.- n-Decane
%Imput data
R = 83.1439;%cm^3-bar/K-gmol
P = 6.8;%bar
T = 71.1;%C
T = T+273.15;%K
zi = [.35 .45 .2];
Tc = [190.6 425.2 617.7];%K
Tr = T./Tc;
Pc = [4.54 3.8 2.12];%Mpa
Pc = Pc*10;%bar
wi = [.008 .199 .489];
Mi = [16.04 58.12 142.29];%g/mole
MT = 100;%kg
MT = MT*1000;%g

E = [1e-6 1e-6 1e-6];
dKi = 1000;
Ki = (Pc./P).*exp(5.37.*(1+wi).*(1-(Tc./T)));

while(abs(dKi)>=E);
%Rachford-Rice Equation 
syms x;
eqn = sum((zi.*(1-Ki))./(Ki+(1-Ki).*x)) == 0;
solx = vpasolve(eqn,x,[0,1]);
l = double(solx);
xi = zi./(l+Ki.*(1-l));
yi = Ki.*zi./(l+Ki.*(1-l));

% MODIFIED PENG-ROBINSON (PR)  EOS
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
aV(i) = sum(yi(i).*yi.*((ai(i).*ai).^.5).*(1-kij(i,:)));
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
solnPRV = roots(EqV)
ZV = solnPRV(1);
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


if(fugLi~=fugVi)
Ki = (fugLi./fugVi).*Ki;
dKi = (fugLi./fugVi)-1;
end

end
v = 1-l;

%Results
fprintf('\n');
fprintf('**********************************');
fprintf('\n');
fprintf('      MODIFIED PENG-ROBINSON EOS');
fprintf('\n');
fprintf('\n');
fprintf(' LIQUID ');
fprintf('\n');
fprintf('Flash "l" = %d',l);
fprintf('\n');
fprintf(' Compositions: ');
fprintf('\n');
fprintf('x_methane = %d',xi(1));
fprintf('\n');
fprintf('x_n-butane = %d',xi(2));
fprintf('\n');
fprintf('x_n-decane = %d',xi(3));
fprintf('\n');
fprintf('\n');
fprintf('Molar Volume = %d',VL);
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf(' VAPOR ');
fprintf('\n');
fprintf('Flash "v" = %d',v);
fprintf('\n');
fprintf('\n');
fprintf(' Compositions: ');
fprintf('\n');
fprintf('y_methane = %d',yi(1));
fprintf('\n');
fprintf('y_n-butane = %d',yi(2));
fprintf('\n');
fprintf('y_n-decane = %d',yi(3));
fprintf('\n');
fprintf('\n');
fprintf('Molar Volume = %d',VV);
fprintf('\n');
fprintf('\n');
fprintf('\n');
fprintf('**********************************');
fprintf('\n');
fprintf('\n');


end

