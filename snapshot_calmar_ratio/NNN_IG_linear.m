function [xk, ksi] = NNN_IG_linear(weeks, wk_return_d1, rf, miu) 
%weeks(>N/2),wk_return_d1: whole data，rb: whole data，miu: mean
wk_return_train = wk_return_d1(:, 1:weeks-1); %online
[M, ~] = size(wk_return_train); 

%改成线性的 linprog parameter
f = [-1*miu', rf(weeks), zeros(1, weeks-1)]; % 分别对应y ita ksi(训练数据weeks-1周)
A1 = diag(ones(1,weeks),0)+diag(-1*ones(1,weeks-1),1);   % 依然是weeks个方程
A1 = [zeros(weeks, M+1),A1];                              % 总共weeks个方程（因为ksi代表每个周的dd）
A1 = A1(1:end-2, 1:end-1);                                       % 去掉最后两行 因为第m周不能再减了
A2 = zeros(weeks-1, M+1+weeks-1);
for j = 1:weeks-1
    A2(j, :) = [-sum(wk_return_d1(:, 1:j),2)', 0 , zeros(1, j-1), 1, zeros(1, weeks-1-j)];
end
A = [-A2;  % ξj−sum(ri'y), j = 1, 2, ..., m; weeks
     A2;  % ξj−sum(ri'y), j = 1, 2, ..., m; weeks
     A1]; %  ksi_j≥ksi_j−1,j=1,2,...,m; weeks-1个约束                                
b = [zeros(weeks-1,1);ones(weeks-1,1);zeros(weeks-2,1)];
Aeq = [ones(1,M),-1,zeros(1,weeks-1)];
beq = 0;
lb = zeros(M+1+weeks-1,1);              % 线性不等式约束

%linprog
op = optimoptions('linprog','display','none');
result = linprog(f, A, b, Aeq, beq, lb, [], op);                          % 求解线性规划问题

yk_32 = result(1:M, :);
ita_32 = result(M+1, :);
ksi_32 = result(M+2:end, :);
xk_32 = yk_32/ita_32;

xk = xk_32;
ksi = ksi_32;
end

