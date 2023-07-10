function [CDD_VaR0] = CDD_Var_p_cal(weeks, wk_return_d1, xt_all, theta)
%wk_return_d1:the whole dataset, rb:vector

MyValue = 1000;
options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue,...
    'ConstraintTolerance',1e-6,'MaxIterations',800,'StepTolerance',1e-8);
[~, CDD_VaR0, a] = fmincon(@(rho) CDD_Var_p(weeks, wk_return_d1, xt_all, rho, theta), 0.1, [], [], [], [], 0, [], [], options);

%[~, CDD_VaR0, a] = fmincon(@(rho) EVaR_p(weeks, wk_return_d1, xt_all, rho, theta), 0.01, [], [], [], [], 0, [], []);

if(a ~= 1)
    disp(['p_cal:' num2str(a)])
end
end 