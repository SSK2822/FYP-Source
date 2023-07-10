function [xk, rho, ksi, nolineark] = NNN_IG(weeks, wk_return_d1, rf, miu) 
%weeks(>N/2),wk_return_d1: whole data，rb: whole data，miu: mean
wk_return_train = wk_return_d1(:, 1:weeks-1); %online
[M, ~] = size(wk_return_train); 

% 验证sortino ratio公式34
x0 = 1 / M * ones(M, 1); %initial as equally weighted
y0 = x0 / (x0'*cov(wk_return_train')*x0)^0.5;
ita0 = 1 / (x0'*cov(wk_return_train')*x0)^0.5; 
theta = 0.95;
rho0 = 0.1;
ksi0 = 0.01 * ones(weeks,1);

Aeq = [ones(1,M),-1,0,zeros(1,weeks)]; % 分别对应y ita rho ksi
beq = 0;
A1 = diag(ones(1,weeks),0)+diag(-1*ones(1,weeks-1),1);   % 依然是weeks个方程
A = [zeros(weeks, M+2),A1];                              % 总共weeks个方程（因为ksi代表每个周的dd）
A = A(1:end-1, :);                                       % 去掉最后一行 因为第m周不能再减了
b = zeros(weeks-1, 1);
lb = zeros(M+1+1+weeks,1);              % 线性不等式约束
nonlcon = @(u)nonlinear_constrain(u, weeks, wk_return_d1, theta, rf);         % 非线性二次约束
% [c,ceq] = nonlinear_constrain([y0;ita0;rho0], weeks-1, wk_return_d1, theta, rf)

%fmincon搜索 patternsearch搜索（一般是非凸）
MyValue = 3000;
options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue,...
    'ConstraintTolerance',1e-6,'MaxIterations',1000,'StepTolerance',1e-10);
[result, ~, a] = fmincon(@(u) odepseudo33(u, miu, rf, weeks, wk_return_train), [y0;ita0;rho0;ksi0], A, b, Aeq, beq, lb, [], nonlcon, options);
%[result, ~, a] = patternsearch(@(u) odepseudo33(u, miu, rf, weeks, wk_return_train), [y0;ita0], [], [], A, b, lb, [], nonlcon, options);

%{
%ga 遗传算法 搜索
%options=gaoptimset('Generations',N/2);
%[result, ~, a] = ga(@(u) odepseudo33(u, miu, rf, weeks, wk_return_train),M+1, [], [], A, b, lb, [], nonlcon, options);

%globalsearch搜索 multistart搜索
%obj = @(u)odepseudo33(u, miu, rf, weeks, wk_return_train);
%problem = createOptimProblem('fmincon','x0',[y0;ita0],'objective',obj,'Aineq',[],'bineq', [],...
%     'Aeq',A, 'beq', b,'lb', lb,'ub',[],'nonlcon', nonlcon, 'options', options);
%[result,~] = run(MultiStart,problem,10);
%[result,~] = run(GlobalSearch,problem);

%}


yk_32 = result(1:M, :);
ita_32 = result(M+1, :);
rho_32 = result(M+2, :);
ksi_32 = result(M+3:end, :);
xk_32 = yk_32/ita_32;

xk = xk_32;
rho = rho_32;
ksi = ksi_32;
[nolineark, ~] = nonlinear_constrain(result, weeks, wk_return_d1, theta, rf);


end

