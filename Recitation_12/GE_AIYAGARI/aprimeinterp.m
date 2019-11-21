function aprimeval = aprimeinterp(aold,aprime0,aval,ect)

%interpolates the mapping aold(:,ect) --> aprime0 at the value aval,
%choosing the minimum and maximum asset values when out of bounds

anum = length(aprime0);

if (aval>=aold(1,ect))&&(aval<=aold(anum,ect))
    aprimeval = interp1q(aold(:,ect),aprime0,aval);
elseif (aval<aold(1,ect))
    aprimeval = aprime0(1);
elseif (aval>aold(anum,ect))     
    aprimeval = aprime0(anum);
end
 
end
                