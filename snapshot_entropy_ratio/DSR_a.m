function [DSR0] = DSR_a(weeks, wk_return_d1, x, rb)
%wk_return_d1:the whole dataset, rb:vector
%DSR used for ex-ante version
DSR1 = zeros(1, weeks-1);

for week = 1:weeks-1 
    rj = wk_return_d1(:, week); 
    if (rj'*x>=rb(week))
       DSR1(week)=0;
    elseif (rj'*x<rb(week))   
        DSR1(week) = (rj'*x-rb(week))^2;
    end 
    
end 
DSR0 = sum(DSR1, 2)/(weeks-1);


end 

