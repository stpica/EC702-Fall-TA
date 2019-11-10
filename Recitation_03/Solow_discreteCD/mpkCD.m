function mpk = mpkCD(K, L, A, alpha)
% This function recalls the expression for the marginal product of capital
% for a Coob Douglas production function

mpk = A * alpha * K ^ (alpha - 1) * L ^ (1 - alpha);

end