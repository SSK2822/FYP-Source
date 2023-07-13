function [CDD_VaR0] = CDD_Var_a_cal(weeks, wk_return_d1, x, theta)
MyValue = 500;
options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue,...
    'ConstraintTolerance',1e-6,'MaxIterations',400,'StepTolerance',1e-8);

[~, CDD_VaR0, a] = fmincon(@(rho)CDD_Var_a(weeks, wk_return_d1, x, rho, theta), ...
    0.01, [], [], [], [], 0, [], [], options...
            );
if(a ~= 1)
    disp(['a_cal:' num2str(a)]);
end
end 