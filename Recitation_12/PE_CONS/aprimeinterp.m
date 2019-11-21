function aprimeval = aprimeinterp(aold,aprime0,aval,yct)

%interpolates the mapping aold(:,yct) --> aprime0 at the value aval,
%choosing the minimum and maximum asset values when out of bounds

anum = length(aprime0);

if (aval>=aold(1,yct))&&(aval<=aold(anum,yct))
    aprimeval = interp1q(aold(:,yct),aprime0,aval);
elseif (aval<aold(1,yct))
    aprimeval = aprime0(1);
elseif (aval>aold(anum,yct))     
    aprimeval = aprime0(anum);
end
 
end
                