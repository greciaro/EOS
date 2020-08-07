
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
solnPRV = roots(EqV);
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


%Results
fprintf('\n');
fprintf('**********************************');
fprintf('\n');
fprintf('      PENG-ROBINSON EOS');
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

