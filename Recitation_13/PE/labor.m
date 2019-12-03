function nvalue = labor(zval,kval,alpha,nu,W)
    nvalue = (zval^(1/(1-nu))) * ( (nu/W)^(1/(1-nu)) ) * (kval^(alpha/(1-nu)));
end