function [c,ceq] = nonlinear_constrain(u, weeks, wk_return_d1, theta, rf)
%wk_return_d1:the whole dataset, rb:vector
%DSR used for ex-ante version
[M, ~] = size(wk_return_d1);
y = u(1:M,:);
ita = u(M+1,:);
rho = u(M+2,:);
exp_sum = 0;

for week = 1:weeks-1
    rj = wk_return_d1(:, week); 
    exp_sum = exp_sum + exp(-rj'*y/rho);
end
c = rho * log(1/(weeks-1)*exp_sum) - rho * log(1-theta) - 1;
A = [ones(1,M),-1,0];             % 线性等式约束
ceq = A*u;                        %线性约束
end 
