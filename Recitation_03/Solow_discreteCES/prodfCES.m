function y = prodfCES(K, L, A, sigma, alpha_K, alpha_L)
% This function recalls the CES production function

y = A * (alpha_K ^ (1 / sigma) * K ^ (1 - 1 / sigma) + alpha_L ^ (1 / sigma) * L ^ (1 - 1 / sigma)) ^ (sigma / (sigma - 1));
    
end