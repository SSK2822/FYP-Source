function obj = odepseudoObjective(u, wk_return_train)
[M, ~] = size(wk_return_train);
y = u(1:M, :);
ita = u(M+1, :);
obj = (y'*cov(wk_return_train')*y)^0.5;
end