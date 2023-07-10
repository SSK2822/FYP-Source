function obj = odepseudo33(u, miu, rf, weeks, wk_return_train)
%fmincon求解器

[M, ~] = size(wk_return_train);
y = u(1:M, :);
ita = u(M+1, :);
obj = -(miu'*y-rf(weeks)*ita); % 注意公式33原式是max，所以要加负号求最小值


end