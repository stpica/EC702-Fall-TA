function mpk=mpkCES(K, L, A, sigma, alpha_K, alpha_L)
% This function recalls the expression for the marginal product of capital

mpk=A*alpha_K^(1/sigma)*K^(-1/sigma)*(alpha_K^(1/sigma)*K.^(1-1/sigma)+alpha_L^(1/sigma)*L.^(1-1/sigma)).^(1/(sigma-1));
    
end