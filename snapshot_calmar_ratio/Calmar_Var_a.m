function [Max_Drawdown0_a] = Calmar_Var_a(weeks, wk_return_d1, x)
%wk_return_d1:the whole dataset, rb:vector
%DSR used for ex-ante version
Drawdown = zeros(1, weeks-1); 

for week = 1:weeks-1
    rj = wk_return_d1(:, 1:week); 
    accu_return = zeros(1, week);
    for i = 1:week
        rkk = wk_return_d1(:, 1:i);
        accu_return(i) = sum(rkk'*x) ;
    end
    Drawdown(week) = max(accu_return)-sum(rj'*x); %week期的maximum drawndown
end
Max_Drawdown0_a = max(Drawdown);

end 