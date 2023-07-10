function [xk, rho, nolineark] = NNN_IG(weeks, wk_return_d1, rf, miu) 
%weeks(>N/2),wk_return_d1: whole data，rb: whole data，miu: mean
wk_return_train = wk_return_d1(:, 1:weeks-1); %online
[M, ~] = size(wk_return_train); 

% 验证sortino ratio公式34
x0 = 1 / M * ones(M, 1); %initial as equally weighted
y0 = x0 / (x0'*cov(wk_return_train')*x0)^0.5;
ita0 = 1 / (x0'*cov(wk_return_train')*x0)^0.5; 
theta = 0.95;
rho0 = 0.1;
lb = zeros(M+2,1);              % 线性不等式约束
nonlcon = @(u)nonlinear_constrain(u, weeks, wk_return_d1, theta, rf);         % 非线性二次约束
% [c,ceq] = nonlinear_constrain([y0;ita0;rho0], weeks-1, wk_return_d1, theta, rf)

%fmincon搜索 patternsearch搜索（一般是非凸）
MyValue = 10e6;
options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue,...
    'ConstraintTolerance',1e-7,'MaxIterations',5000,'StepTolerance',1e-14);
[result, ~, a] = fmincon(@(u) odepseudo33(u, miu, rf, weeks, wk_return_train), [y0;ita0;rho0], [], [], [], [], lb, [], nonlcon, options);
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
rho_32 = result(M+2,:);
xk_32 = yk_32/ita_32;

xk = xk_32;
rho = rho_32;
[nolineark, ~] = nonlinear_constrain(result, weeks, wk_return_d1, theta, rf);


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
%ode hyperparameter 常微分方程参数
eps1 = 0.0001;%epsilon_u
eps2 = 0.0001;
alpha = 2;
beta = 2;
t_end = 0.3;
tspan = [0 t_end];
%tspan = 0:step:t_end;

%ode solver(sometimes use ode15s if needed)
%[t,solution] = ode45(@(t,u) odepseudoNN(u, P1, q1, I1, miu, wk_return_train, ep1), tspan, [x0;ksi0]);
[t,solution] = ode15s(@(t,u) odepseudoONENN(u, miu, wk_return_train, rf, weeks, eps1, eps2, alpha, beta), tspan, [y0;ita0]);
solution0 = solution(end, :)'; %取最后平稳的那个solution
y0 = solution0(1:M);
ita0 = solution0(M+1);
xk_onenn = solution0(1:M)/solution0(M+1)
equ1 = sum(xk_onenn)
equ2 = (miu'*y0) - rf(weeks)*ita0
[m1, ~] = size(solution);
sr_plot = zeros(1, m1);
ita_plot = zeros(1, m1);
for j = 1:m1
    u_plot = solution(j, :)'; % 注意是M+1维
    y_plot = u_plot(1:M);
    ita_plot(j) = u_plot(M+1);
    x_plot = y_plot/ita_plot(j); % y/ita变换得到x
    ratio_a_my = 1/((y_plot'*cov(wk_return_train')*y_plot)^0.5); %as defined in 34 应该直接用这个式子 sp
    sr_plot(j) = ratio_a_my;
end
sr_plot(end)

figure(1)
p1 = plot(t, sr_plot, '-','linewidth', 1.5, 'markersize', 5,'color', 'r');
hold on

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
for j = 1:m
    u_plot_pseudo = u_pseudo(j, :);
    x_plot_pseudo = max(0, u_plot_pseudo);
    x_plot_pseudo = x_plot_pseudo';
    My_wk_rt_a_temp = x_plot_pseudo' * mean(wk_return_train, 2) - rf(weeks);
    ratio_a_my = My_wk_rt_a_temp / (x_plot_pseudo'*cov(wk_return_train')*x_plot_pseudo)^0.5; %as defined in sr_a
    sr_plot_pseudo(j) = ratio_a_my;
end
t_plot_pseudo = t_pseudo;
linewidth_gca = 1.6;
fontsize_gca = 24;
fontsize_label = 28;
set(gca, 'linewidth', linewidth_gca, 'fontsize', fontsize_gca);
ylabel('\itSR', 'fontsize', fontsize_label)   ;
xlim([-0.002 0.1]);
ylim([0, 0.2])
title('FTSE'); 
hold on
p2 = plot(t_plot_pseudo, sr_plot_pseudo, '-.','linewidth', 1.5, 'markersize', 5,'color', 'b'); %直接算sharp ratio
set(gcf, 'unit', 'centimeters', 'position', [10, 10, 16, 10])
fontsize_legend = 20;
p1 = [];
h = legend(p1, ['One step']); % ,sprintf('\n'), 'convex optimization'
legend('boxoff')
ah=axes('position',get(gca,'position'),'visible','off');
h2 = legend(ah, p2, 'Pseudoconvex'); %  Optimization
legend('boxoff')
set(h,'fontsize',fontsize_legend);
set(h2,'fontsize',fontsize_legend);
hold on

sr_plot_pseudo(end)
xk = xk_pseudo %output

%}
end

