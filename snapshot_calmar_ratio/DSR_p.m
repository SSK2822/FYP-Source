function [DSR0] = DSR_p(weeks, wk_return_d1, x, rb)
%whole dataset input, x is decision matrix(multi-period)
%used for ex-post version
[~, N] = size(wk_return_d1);
wk_return = wk_return_d1(:, N/2+1:end);
%wk_return = wk_return_d1(:, N/2:end);
DSR1 = zeros(1, weeks-N/2); 

for week = 1:(weeks-N/2) %m in paper
    rj = wk_return(:, week);
    xj = x(:, week);
    rfj = rb(week+N/2);
    if (rj'*xj>=rfj)
       DSR1(week)=0;
    elseif (rj'*xj<rfj)  
        DSR1(week) = (rj'*xj-rfj)^2;
    end 
    
end 
DSR0 = sum(DSR1, 2)/(weeks-N/2);




end %输出是一个值

