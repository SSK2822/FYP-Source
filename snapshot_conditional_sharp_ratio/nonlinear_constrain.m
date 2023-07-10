function [c,ceq] = nonlinear_constrain(u, weeks, wk_return_train, rf, theta)
[M, ~] = size(wk_return_train);
y = u(1:M, :);
ita = u(M+1, :);
rho = u(M+2, :);
sig = zeros(1, weeks-1);
rb = rf(weeks); % 注意目标函数不能用rb 每周都不一样

for week = 1:weeks-1
    rj = wk_return_train(:,week);
    if(-rj'*y-rho>0)
        sig(week) = -rj'*y-rho;
    elseif(-rj'*y-rho<=0)
        sig(week) = 0;
    end
end
A = [ones(1,M),-1,0]; 
c = rho + 1/((weeks-1)*(1-theta))*sum(sig) -1 ;
ceq = A*u;

end