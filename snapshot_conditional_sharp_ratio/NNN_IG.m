function [xk, rhok, nolineark] = NNN_IG(weeks, wk_return_d1, rf, miu) 
%weeks(>N/2),wk_return_d1: whole data，rb: whole data，miu: mean
wk_return_train = wk_return_d1(:, 1:weeks-1); %online
[M, ~] = size(wk_return_train); 
% PNN initial variables
%{
A = [ones(1,M),-1;miu',-rf(weeks)];
b = [0;1];
P1 = (A')*inv(A*A')*A;    % 直接计算新公式 此时P公式的A需要根据约束条件重写
I1 = eye(M+1);
q1 = (A')*inv(A*A')*b;
P2 = 1 / M * ones(M, M);  % 直接计算sortino ratio 的伪凸
I2 = eye(M); 
q2 = 1 / M * ones(M, 1);
%}

% 验证omega ratio公式33
x0 = 1 / M * ones(M, 1); %initial as equally weighted
y0 = x0 / (x0'*cov(wk_return_train')*x0)^0.5;
ita0 = 1 / (x0'*cov(wk_return_train')*x0)^0.5; 
rho0 = 0.01;
theta = 0.95;

lb = zeros(M+2,1);              % 线性不等式约束
% V = cov(wk_return_train');      % 56*56
% v0 = [[V;zeros(1,length(V))],zeros(length(V)+1,1)]; % 二次约束少了一维，要补上
nonlcon = @(u)nonlinear_constrain(u, weeks, wk_return_train, rf, theta);          % 非线性二次约束
obj = @(u)odepseudo33(u, miu, rf, weeks, wk_return_train);


MyValue = 10e6;
options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue,...
    'ConstraintTolerance',1e-8,'MaxIterations',6000,'StepTolerance',1e-14);
[result, ~, a] = fmincon(@(u) odepseudo33(u, miu, rf, weeks, wk_return_train), [x0;1;rho0], [], [], [], [], lb, [], nonlcon, options);

% problem = createOptimProblem('fmincon','x0',[y0;ita0],'objective',obj,'Aineq',[],'bineq', [],...
%      'Aeq',A, 'beq', b,'lb', lb,'ub',[],'nonlcon', nonlcon, 'options', options);
% [result,~] = run(GlobalSearch,problem);

yk_33 = result(1:M, :);
ita_33 = result(M+1, :);
rho_33 = result(M+2, :);
xk_33 = yk_33/ita_33;

xk = xk_33;
rhok = rho_33;
[nolineark, ~] = nonlinear_constrain(result, weeks, wk_return_train, rf, theta);

end

