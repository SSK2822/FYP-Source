function [Martin_Drawdown0_a] = Martin_Var_a(weeks, wk_return_d1, x)
Drawdown = zeros(1, weeks-1); 

for week = 1:weeks-1
    rj = wk_return_d1(:, 1:week); 
    accu_return = zeros(1, week);
    for i = 1:week
        rkk = wk_return_d1(:, 1:i);
        accu_return(i) = sum(rkk'*x) ;
    end
    Drawdown(week) = max(accu_return)-sum(rj'*x); 
    Drawdown(week) = Drawdown(week)^2;
end
Martin_Drawdown0_a = sum(Drawdown)/(weeks-1);

end 