function [EVaR0] = Calmar_Var_p_cal(weeks, wk_return_d1, xt_all, theta)
%wk_return_d1:the whole dataset, rb:vector

MyValue = 10e6;
options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue,...
    'ConstraintTolerance',1e-7,'MaxIterations',5000,'StepTolerance',1e-14);
[~, EVaR0, a] = fmincon(@(rho) EVaR_p(weeks, wk_return_d1, xt_all, rho, theta), 0.1, [], [], [], [], 0, [], [], options);

%[~, EVaR0, a] = fmincon(@(rho) EVaR_p(weeks, wk_return_d1, xt_all, rho, theta), 0.01, [], [], [], [], 0, [], []);

if(a ~= 1)
    disp(['p_cal:' num2str(a)])
end
end 