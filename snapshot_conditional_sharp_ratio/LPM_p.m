function [LPM0] = LPM_p(weeks, wk_return_d1, x, rb)
%whole dataset input, x is decision matrix(multi-period)
%used for ex-post version
[~, N] = size(wk_return_d1);
wk_return = wk_return_d1(:, N/2+1:end); % 后半部分数据
%wk_return = wk_return_d1(:, N/2:end);
LPM1 = zeros(1, weeks-N/2); 

for week = 1:(weeks-N/2) %m in paper
    rj = wk_return(:, week);
    xj = x(:, week);
    rfj = rb(week+N/2); %写的不对吧
    if (rj'*xj>=rfj)
       LPM1(week)=0;
    elseif (rj'*xj<rfj)  
        LPM1(week) = (rfj-rj'*xj);
    end 
    
end 
LPM0 = sum(LPM1, 2)/(weeks-N/2);




end %输出是一个值

