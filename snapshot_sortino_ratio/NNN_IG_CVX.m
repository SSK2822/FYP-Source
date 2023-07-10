function [xk, nolineark] = NNN_IG_CVX(weeks, wk_return_d1, rf, miu) 
%weeks(>N/2),wk_return_d1: whole data，rb: whole data，miu: mean
wk_return_train = wk_return_d1(:, 1:weeks-1); %online
[M, ~] = size(wk_return_train); 

x0 = 1 / M * ones(M, 1); %initial as equally weighted
y0 = x0 / (x0'*cov(wk_return_train')*x0)^0.5;
ita0 = 1 / (x0'*cov(wk_return_train')*x0)^0.5; 
% nonlinear_1 = nonlinear_constrain([y0;ita0], weeks, wk_return_train, rf)

A = [ones(1,M),-1];             % 线�?�等式约�?
b = 0;
lb = zeros(M+1,1);              % 线�?�不等式约束
V = cov(wk_return_train');      % 56*56
v0 = [[V;zeros(1,length(V))],zeros(length(V)+1,1)]; % 二次约束少了�?维，要补�?
nonlcon = @(u)nonlinear_constrain(u, weeks, wk_return_train, rf);          % 非线性二次约�?

%fmincon搜索 patternsearch搜索（一般是非凸�?
%MyValue = 10e6;
%options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue,...
%    'ConstraintTolerance',1e-6,'MaxIterations',5000,'StepTolerance',1e-12);
%[result, ~, a] = fmincon(@(u) odepseudo33(u, miu, rf, weeks, wk_return_train), [x0;1], [], [], [], [], lb, [], nonlcon, options);
%[result, ~, a] = patternsearch(@(u) odepseudo33(u, miu, rf, weeks, wk_return_train), [y0;ita0], [], [], A, b, lb, [], nonlcon, options);

% variable
u = sdpvar(M+1,1,'full');%uΪʵ��
y = u(1:M, :);
ita = u(M+1, :);

% constrain
A = [ones(1,M),-1]; 
rb = rf(weeks);
wk_return_train_a = wk_return_train(:,1:weeks-1);
c = [];
c = [c;A*u==0;u>=0];
c = [c;sum(max((rb*ita-wk_return_train_a'*y),0).^2)/(weeks-1) <= 1];
% c = [c;(sum(dsr,2)/(weeks-1))<=1];

%{
for week = 1:weeks-1
    rj = wk_return_train(:,week);
    cdsr(week) = (max(rb*ita-rj'*y,0).^2);
%    d = binvar(2,1);
%    c = [c, sum(d) == 1,...
%        implies(d(1), [rj'*y>=rb*ita, dsr(week) == 0]);...
%        implies(d(2), [rj'*y<=rb*ita, dsr(week) == (rj'*y-rb*ita)^2])];
end
%}


% objective function
obj = -(miu'*y-rf(weeks)*ita);
% select solver
ops = sdpsettings('verbose',2,'solver','cplex'); % 'maxtime',5555
ops.bnb.maxtime=20;
result = optimize(c,obj,ops);%��� 
if result.problem == 0
    result = value(u)
    obj_result = value(obj)
else
    disp('������')
end 

size(result)
yk_32 = result(1:M, :);
ita_32 = result(M+1, :);
xk_32 = yk_32/ita_32;

xk = xk_32;

[nolineark, ~] = nonlinear_constrain(result, weeks, wk_return_train, rf);
end

