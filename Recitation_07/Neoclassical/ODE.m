%input of system of ordinary differential equations (as a vector)
 
funcODE=[ A_neutral*k^alpha*(ell)^(1-alpha)-c-(nPop+delta+gx)*k; ...         %ODE: k
          c*(theta*((1-mu)*alpha*A_neutral*k^(alpha-1)*ell^(1-alpha)-rho-delta) -gx);  %ODE: c
             ];
         
         %first state, then control