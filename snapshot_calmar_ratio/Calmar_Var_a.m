function [Max_Drawdown0_a] = Calmar_Var_a(weeks, assets_retur, x)
Drawdown = zeros(1, weeks-1); 

for week = 1:weeks-1
    rj = assets_retur(:, 1:week); 
    accu_return = zeros(1, week);
    for i = 1:week
        rkk = assets_retur(:, 1:i);
        accu_return(i) = sum(rkk'*x) ;
    end
    Drawdown(week) = max(accu_return)-sum(rj'*x); 
end
Max_Drawdown0_a = max(Drawdown);

end 