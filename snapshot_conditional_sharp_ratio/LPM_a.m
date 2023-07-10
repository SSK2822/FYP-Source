function [LPM0] = LPM_a(weeks, wk_return_d1, x, rb)
%wk_return_d1:the whole dataset, rb:vector all
%DSR used for ex-ante version
LPM1 = zeros(1, weeks-1);

for week = 1:weeks-1 
    rj = wk_return_d1(:, week); 
    if (rj'*x>=rb(week))
       LPM1(week)=0;
    elseif (rj'*x<rb(week))   
        LPM1(week) = (rb(week)-rj'*x);
    end 
    
end 
LPM0 = sum(LPM1, 2)/(weeks-1);


end 

