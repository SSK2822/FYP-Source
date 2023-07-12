function [portfolio_weights, ksi] = NNN_IG_linear(week, asset_returns, risk_free_rate, mean_return) 
asset_returns_train = asset_returns(:, 1:week-1); % online
[num_assets, ~] = size(asset_returns_train); 

% Change to linear linprog parameter
f = [-1*mean_return', risk_free_rate(week), zeros(1, week-1)]; 
A1 = diag(ones(1,week),0)+diag(-1*ones(1,week-1),1);   
A1 = [zeros(week, num_assets+1),A1];                              
A1 = A1(1:end-2, 1:end-1);                                       
A2 = zeros(week-1, num_assets+1+week-1);
for j = 1:week-1
    A2(j, :) = [-sum(asset_returns(:, 1:j),2)', 0 , zeros(1, j-1), 1, zeros(1, week-1-j)];
end
A = [-A2;  % ξj−sum(ri'y), j = 1, 2, ..., m; weeks
     A2;  % ξj−sum(ri'y), j = 1, 2, ..., m; weeks
     A1]; %  ksi_j≥ksi_j−1,j=1,2,...,m; weeks-1 constraints                                
b = [zeros(week-1,1);ones(week-1,1);zeros(week-2,1)];
Aeq = [ones(1,num_assets),-1,zeros(1,week-1)];
beq = 0;
lb = zeros(num_assets+1+week-1,1);              

%linprog
op = optimoptions('linprog','display','none');
result = linprog(f, A, b, Aeq, beq , lb , [], op);                          

yk_32 = result(1:num_assets , :);
ita_32 = result(num_assets+1 , :);
ksi_32 = result(num_assets+2:end , :);
xk_32 = yk_32/ita_32;

portfolio_weights = xk_32;
ksi = ksi_32;
end
