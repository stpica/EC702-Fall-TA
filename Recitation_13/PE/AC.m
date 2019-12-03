function ACval = AC(kval,kprimeval,delta,gamma)
ival = kprimeval - (1-delta)*kval;
ACval = (gamma/2) * ( ival^2);
end