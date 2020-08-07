
% For Vapor
for i=1:length(ai);
aV(i) = sum(yi(i).*yi.*((ai(i).*ai).^(1/2)));
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
solnRKV = roots(EqV);
ZV = solnRKV(1);
%Volume
VV = ZV*R*T/P;

%Fugacity
deltaVi = 2.*(ai./aV).^.5;
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
fprintf('      REDLICH-KOWNG EOS');
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
