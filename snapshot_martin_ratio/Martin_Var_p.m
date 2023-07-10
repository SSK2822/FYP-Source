function [Martin_Drawdown0] = Martin_Var_p(weeks, wk_return_d1, xt_all)
%whole dataset input, x is decision matrix(multi-period)
%used for ex-post version
%weeks: the ith weeks of test data
%计算后验ADD
[~, N] = size(wk_return_d1);
wk_return = wk_return_d1(:, N/2+1:end);
Drawdown = zeros(1, weeks); 

for week = 1:weeks
    rj = wk_return(:, 1:week);
    accu_return = zeros(1, week);
    for i = 1:week
        rkk = wk_return(:, 1:i);
        xtt = xt_all(:, 1:i);
        accu_return(i) =  sum(diag(rkk'*xtt));
    end
    
    Drawdown(week) = max(accu_return)-sum(diag(rj'*xt_all(:, 1:week))); %week期的maximum drawndown 一定是大于0的
    Drawdown(week) = Drawdown(week)^2; %注意Martin ratio要加平方
end
Martin_Drawdown0 = sum(Drawdown)/weeks;

end

