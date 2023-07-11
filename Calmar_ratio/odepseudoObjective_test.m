function obj = odepseudoObjective(u, rf, weeks, wk_return_train)
[M, ~] = size(wk_return_train);
y = u(1:M, :);
ita = u(M+1, :);
dsr = zeros(1, weeks-1);

for week = 1:weeks-1
    rj = wk_return_train(:,week);
    if (rj'*y>=rf(week))
        dsr(week)=0;
    elseif (rj'*y<rf(week))
        dsr(week) = (rj'*y-rf(week))^2;
    end
end

obj = (sum(dsr,2)/(weeks-1))^0.5; %注意公式37是min
end