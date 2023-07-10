function [c,ceq] = nonlinear_constrain(u, weeks, wk_return_d1, theta, rf)
%wk_return_d1:the whole dataset, rb:vector
%DSR used for ex-ante version
[M, ~] = size(wk_return_d1);
y = u(1:M,:);
ita = u(M+1,:);
ksi = u(M+2:end, :); %ksi总共有weeks个变量

c = zeros(1, 2*(weeks-1));
con_sum = zeros(1, weeks-1);

for week = 1:weeks-1
    rj = wk_return_d1(:, 1:week); 
    con_sum(week) = max(ksi(week)-sum(rj'*y),0);
end
c = sum(con_sum)/(weeks-1)-1;
ceq = [];

end 
