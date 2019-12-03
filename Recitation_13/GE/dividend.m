function dval = dividend(zval,kval,kprimeval,alpha,nu,W,delta,gamma)
 
nval = labor(zval,kval,alpha,nu,W);
yval = output(zval,kval,nval,alpha,nu);
ival = kprimeval - (1-delta)*kval;
ACval = AC(kval,kprimeval,delta,gamma);
dval = (1-nu)*yval - ival - ACval;



end