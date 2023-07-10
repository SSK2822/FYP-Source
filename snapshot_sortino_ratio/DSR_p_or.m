function [DSR0] = DSR_p_or(weeks, wk_return_d1, x, rb)%r是矩阵，x是列向量，rb是行向量
%后验DSR 训练用%weeks是测试集的第m周
[M, N] = size(wk_return_d1);
%wk_return = wk_return_d1(:, N-468:end);
wk_return = wk_return_d1(:, N/2:end);
DSR1 = zeros(1, weeks); %初始化一个中间行向量用来求均值

for week = 1:weeks %weeks是公式里的第m周
    rj = wk_return(:, week);%rj是列向量
    if (rj'*x>=rb(week+N/2))
    %if (rj'*x>=rb(week+N-469))
       DSR1(week)=0;
    %elseif (rj'*x<rb(week+N-469))  
    elseif (rj'*x<rb(week+N/2))  
        %DSR1(week) = (rj'*x-rb(week+N-469))^2;
        DSR1(week) = (rj'*x-rb(week+N/2))^2;
    end %循环2
    
end %循环1
DSR0 = sum(DSR1, 2)/weeks;




end %输出是一个值

