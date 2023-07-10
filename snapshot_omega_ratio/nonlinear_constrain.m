function [c,ceq] = nonlinear_constrain(u, weeks, wk_return_train, rf)
[M, ~] = size(wk_return_train);
y = u(1:M, :);
ita = u(M+1, :);
lpm = zeros(1, weeks-1);
rb = rf(weeks); % 注意目标函数不能用rb 每周都不一样

for week = 1:weeks-1
    rj = wk_return_train(:,week);
    if (rj'*y>=rb*ita)
        lpm(week)=0;
    elseif (rj'*y<rb*ita)
        lpm(week) = rb*ita-rj'*y;
    end
end

A = [ones(1,M),-1]; 
c = (sum(lpm,2)/(weeks-1)) - 1; %公式33非线性约束
ceq = A*u;

end