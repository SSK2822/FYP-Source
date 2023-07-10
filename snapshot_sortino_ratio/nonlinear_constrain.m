function [c,ceq] = nonlinear_constrain(u, weeks, wk_return_train, rf)
[M, ~] = size(wk_return_train);
y = u(1:M,:);
ita = u(M+1,:);
dsr = zeros(1, weeks-1);
rb = rf(weeks);

for week = 1:weeks-1
    rj = wk_return_train(:,week);
    if (rj'*y>rb*ita)
        dsr(week)=0;
    elseif (rj'*y<=rb*ita)
        dsr(week) = (rj'*y-rb*ita)^2;
    end
end

A = [ones(1,M),-1]; 
c = (sum(dsr,2)/(weeks-1))-1; %公式32非线性约束 注意不用加根号
ceq = A*u;                    %线性约束
end