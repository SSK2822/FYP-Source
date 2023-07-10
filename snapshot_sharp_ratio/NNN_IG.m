function [xk, ratiok, nolineark] = NNN_IG(weeks, wk_return_d1, rf, miu) 
%weeks(>N/2),wk_return_d1: whole data，rb: whole data，miu: mean
wk_return_train = wk_return_d1(:, 1:weeks-1); %online
[M, ~] = size(wk_return_train); 
A = [ones(1,M),-1;miu',-rf(weeks)];
b = [0;1];
P1 = (A')*inv(A*A')*A;    % 直接计算新公式 此时P公式的A需要根据约束条件重写
I1 = eye(M+1);
q1 = (A')*inv(A*A')*b;
P2 = 1 / M * ones(M, M); % 直接计算sharp ratio 的伪凸
I2 = eye(M); 
q2 = 1 / M * ones(M, 1);

%{
% 验证公式35
x0 = 1 / M * ones(M, 1); %initial as equally weighted
y0 = x0 / (miu'*x0-rf(weeks));
ita0 = 1 / (miu'*x0-rf(weeks)); %FTSE

MyValue = 10e4;
options = optimoptions('fmincon', 'Algorithm','sqp', 'MaxFunctionEvaluations',MyValue);
[result, ~, a] = fmincon(@(u) odepseudoObjective(u, wk_return_train), [y0;ita0], [], [], A, b, zeros(M+1, 1), [], [], options);
yk_convex = result(1:M, :);
ita_convex = result(M+1, :);
xk_convex = yk_convex/ita_convex
ratio_a_my_convex1 = 1/((yk_convex'*cov(wk_return_train')*yk_convex)^0.5) %as defined in 34 应该直接用这个式子 sp
%}

% 验证公式33
x0 = 1 / M * ones(M, 1); %initial as equally weighted
y0 = x0 / (x0'*cov(wk_return_train')*x0)^0.5;
ita0 = 1 / (x0'*cov(wk_return_train')*x0)^0.5; 

A = [ones(1,M),-1];             % 线性等式约束
b = 0;
lb = zeros(M+1,1);              % 线性不等式约束
V = cov(wk_return_train');   % 56*56
v0 = [[V;zeros(1,length(V))],zeros(length(V)+1,1)]; % 二次约束少了一维，要补上
nonlcon = @(u)nonlinear_constrain(u,v0);            % 非线性二次约束

MyValue = 500000;
options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue);
[result, ~, a] = fmincon(@(u) odepseudo33(u, miu, rf, weeks, wk_return_train), [y0;ita0], [], [], A, b, lb, [], nonlcon, options);
yk_31 = result(1:M, :);
ita_31 = result(M+1, :);
xk_31 = yk_31/ita_31; % 公式33最终的x解
ratio_a_my_31 = miu'*yk_31 - rf(weeks)*ita_31;

xk = xk_31;
ratiok = ratio_a_my_31;
nolineark = yk_31'*V*yk_31;

%{
% 验证公式45
x0 = 1 / M * ones(M, 1); %initial as equally weighted
y0 = x0 / (x0'*cov(wk_return_train')*x0)^0.5; % 初始值与公式33保持一致

A = ones(1,M);             % 线性等式约束
b = 1;
lb = zeros(M,1);              % 线性不等式约束
v0 = cov(wk_return_train');   % 56*56
nonlcon = @(y)nonlinear_constrain(y,v0);            % 非线性二次约束  

MyValue = 10e4;
options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue);
[result, ~, a] = fmincon(@(y) odepseudo45(y, miu), x0, [], [], A, b, lb, [], nonlcon, options);
yk_45 = result(1:M, :);
ita_45 = 1;
xk_45 = yk_45/ita_45 % 公式33最终的x解
ratio_a_my_45 = miu'*yk_45
%}

%{
% pseudoconvex
step_pseudo = 0.002;
t_end_pseudo = 0.4;
tspan_pseudo = 0:step_pseudo:t_end_pseudo; 
[t_pseudo,u_pseudo] = ode45(@(t,u) odepseudo(u, P2, q2, I2, miu, wk_return_train, 0.0001, rf(weeks)), tspan_pseudo, 1 / M * ones(M, 1));
solution0 = u_pseudo(end, :)'; 
xk_pseudo = max(0, solution0);
[m, ~] = size(u_pseudo);
sr_plot_pseudo = zeros(1, m);

xk = xk_pseudo %output
sr_plot_pseudo(end)
%}

end

