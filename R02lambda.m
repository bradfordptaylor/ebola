function finlambda = R02lambda(R0,rhod,pars)



%number of exposed classes in simulation
boxcont = zeros(1,1+pars.ne);
for jj=0:pars.ne
    boxcont(jj+1)=(pars.ne).^(jj).*(pars.Te).^(pars.ne-jj).*nchoosek(pars.ne,jj);%*(1./pars.Te./pars.ne).^jj
end
rhsmultiply = [pars.Ti.*pars.Td...
    (pars.Td+pars.Ti)...
    1];
rhs1 = [rhsmultiply(1).*boxcont 0 0];
rhs2 = [0 rhsmultiply(2).*boxcont 0];
rhs3 = [0 0 rhsmultiply(3).*boxcont];
lhs =  zeros(size(rhs3));
lhs(end) = R0.*pars.ne^pars.ne;
lhs(end-1) = R0.*pars.ne^pars.ne.*pars.Td.*(1-rhod);
totpoly = rhs1+rhs2+rhs3-lhs;
finlambdatemp = roots(totpoly);
finlambdatemp = finlambdatemp(imag(finlambdatemp)==0);
finlambdatemp = finlambdatemp(finlambdatemp>0);
finlambda = finlambdatemp(1);