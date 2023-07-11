function [c,ceq] = nonlinear_constrain(u, weeks, wk_return_d1, theta, rf)
%wk_return_d1:the whole dataset, rb:vector
%DSR used for ex-ante version
[M, ~] = size(wk_return_d1);
y = u(1:M,:);
ita = u(M+1,:);
ksi = u(M+2:end, :);

c = zeros(1, 2*(weeks-1));
con_sum = zeros(1, weeks-1);

for week = 1:weeks-1
    rj = wk_return_d1(:, 1:week); 
    con_sum(week) = ksi(week)-sum(rj'*y);
    c(week) = con_sum(week)-1;
    c(weeks-1+week) = 0-con_sum(week);
end
% c = [con_sum-1; 0-con_sum];
ceq = [];

end 
