function [CDD_VaR0] = CDD_Var_p_index_cal(weeks, wk_return_b1, theta)
MyValue = 1000;
options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue,...
    'ConstraintTolerance',1e-6,'MaxIterations',800,'StepTolerance',1e-8);
[~, CDD_VaR0, a] = fmincon(@(rho) CDD_Var_p_index(weeks, wk_return_b1, rho, theta), 0.1, [], [], [], [], 0, [], [], options);
if(a ~= 1)
    disp(['p_cal:' num2str(a)])
end
end 