function [portfolio_weights, rho, ksi, nolinear] = NNN_IG(week, asset_returns, risk_free_rate, mean_return) 
asset_returns_train = asset_returns(:, 1:week-1); 
[num_assets, ~] = size(asset_returns_train); 

initial_weights = 1 / num_assets * ones(num_assets, 1); % initial as equally weighted
y0 = initial_weights / (initial_weights'*cov(asset_returns_train')*initial_weights)^0.5;
ita0 = 1 / (initial_weights'*cov(asset_returns_train')*initial_weights)^0.5; 
theta = 0.95;
rho0 = 0.1;
ksi0 = 0.01 * ones(week,1);

Aeq = [ones(1,num_assets),-1,zeros(1,week)];
beq = 0;
A1 = diag(ones(1,week),0)+diag(-1*ones(1,week-1),1);   
A = [zeros(week, num_assets+1),A1];                              
A = A(1:end-1, :);                                       
b = zeros(week-1, 1);
lb = zeros(num_assets+1+week,1);              
nonlcon = @(u)nonlinear_constrain(u, week, asset_returns, theta, risk_free_rate);         
% [c,ceq] = nonlinear_constrain([y0;ita0;rho0], week-1, asset_returns, theta, risk_free_rate)

% fmincon search patternsearch search
MyValue = 4000;
options = optimoptions('fmincon','Algorithm','sqp', 'MaxFunctionEvaluations',MyValue,...
    'ConstraintTolerance',1e-6,'MaxIterations',3000,'StepTolerance',1e-10);
[result, ~, a] = fmincon(@(u) odepseudo33(u, mean_return, risk_free_rate, week, asset_returns_train), [y0;ita0;rho0;ksi0], A, b, Aeq, beq, lb, [], nonlcon, options);
%[result, ~, a] = patternsearch(@(u) odepseudo33(u, mean_return, risk_free_rate, week, asset_returns_train), [y0;ita0], [], [], Aeq, beq , lb , [], nonlcon , options);

%{
%ga Genetic algorithm search
%options=gaoptimset('Generations',num_weeks/2);
%[result , ~ , a ] = ga(@(u) odepseudo33(u , mean_return , risk_free_rate , week , asset_returns_train) , num_assets+1 , [] , [] , Aeq , beq , lb , [] , nonlcon , options);

%globalsearch search multistart search
%obj = @(u)odepseudo33(u , mean_return , risk_free_rate , week , asset_returns_train);
%problem = createOptimProblem('fmincon','x0',[y0;ita0],'objective',obj,'Aineq',[],'bineq', [],...
%     'Aeq',Aeq, 'beq', beq,'lb', lb,'ub',[],'nonlcon', nonlcon, 'options', options);
%[result,~] = run(MultiStart,problem,10);
%[result,~] = run(GlobalSearch,problem);
%}

yk_32 = result(1:num_assets, :);
ita_32 = result(num_assets+1, :);
rho_32 = result(num_assets+2, :);
ksi_32 = result(num_assets+3:end, :);
xk_32 = yk_32/ita_32;
rho = rho_32;
portfolio_weights = xk_32;
ksi = ksi_32;
[nolinear, ~] = nonlinear_constrain(result, week, asset_returns, theta, risk_free_rate);

end
