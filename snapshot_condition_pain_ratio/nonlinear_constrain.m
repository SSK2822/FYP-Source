function [c,ceq] = nonlinear_constrain(u, weeks, wk_return_d1, theta, rf)
%wk_return_d1:the whole dataset, rb:vector
%DSR used for ex-ante version
[M, ~] = size(wk_return_d1);
y = u(1:M,:);
ita = u(M+1,:);
rho = u(M+2,:);
ksi = u(M+3:end, :); %ksi总共有weeks个变量

con_sum = zeros(1, weeks-1);

for week = 1:weeks-1
    rj = wk_return_d1(:, 1:week); 
    con_sum(week) = max(ksi(week)-sum(rj'*y)-rho,0); % 注意Martin ratio约束条件要加平方
end
c = rho + 1/((weeks-1)*(1-theta))*sum(con_sum) - 1; %公式38非线性不等式约束
ceq = [];

end 
