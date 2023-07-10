clc,clear
format long;

filedir1 = './OR_DAX100_CSR/';
filelist1 = dir([filedir1, '*.txt']);
filedir2 = './OR_DAX100_CSR_rho/';
filelist2 = dir([filedir2, '*.txt']);

d1 = load('ORalldata_DAX.txt', ','); %price矩阵(939天数据)
b1 = load('ORb1_DAX.txt', ','); %股指 
rf_92_97 = load('rf_92_97.txt',','); %risk-free

[wk_return_d1, ~] = price2ret(d1', [], 'Periodic'); %simple return
wk_return_d1 = wk_return_d1';
[wk_return_b1, ~] = price2ret(b1, [], 'Periodic'); %simple return

[M, N] = size(wk_return_d1);
rf_test = rf_92_97(N/2+1:N);

My_wk_rt_temp = 0; 
ac_wk_rt_temp = 0;
My_wk_rt = zeros(1, N/2);
ac_wk_rt = zeros(1, N/2);

ew_wk_rt = zeros(1, N/2);
ew_ac_wk_rt = zeros(1, N/2);

x_ew = 1 / M * ones(M, 1);
x_ew_all = 1 / M * ones(M, N/2);
xt_all = zeros(M, N/2);
theta = 0.95;

cvar = zeros(1, N/2);
cvar_a = zeros(1, N/2);
cvar_ew = zeros(1, N/2);
cvar_index = zeros(1, N/2);
csr_a = zeros(1, N/2);
csr = zeros(1, N/2);
csr_ew = zeros(1, N/2);
csr_index = zeros(1, N/2);
index_wk_rt = zeros(1, N/2);
index_ac_wk_rt = zeros(1, N/2);
risk = zeros(1, N/2);
risk_ew = zeros(1, N/2);
risk_index = zeros(1, N/2);
ratio_p_all = zeros(1, N/2);
ew_wk_rt = zeros(1, N/2);
ew_ac_wk_rt = zeros(1, N/2);
ew_wk_rt_temp = 0;
ew_ac_wk_rt_temp = 0;

for i = (N/2+1):N
    wk_return_d1_train = wk_return_d1(:, 1:i-1);%(M, i-1)
    wk_return_d1_test = wk_return_d1(:, i);
    rf = rf_92_97(i);
    miu = mean(wk_return_d1_train, 2);
    filename1 = filelist1(i-N/2);
    filename2 = filelist2(i-N/2);

    x_solution = load(['./OR_DAX100_CSR/',filename1.name]);
    xt_all(:, i-N/2) = x_solution;
    rho_solution = load(['./OR_DAX100_CSR_rho/',filename2.name]);

    %cvar_p
    My_wk_rt_temp = x_solution' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;   
    [rho, cvar0] = CVaR_p_cal(i-N/2, wk_return_d1, xt_all, theta);
    ratio_my = sum(My_wk_rt) / (i-N/2) / cvar0;
    ratio_p_all(i-N/2) = ratio_my;
    

    %cvar_a
    My_wk_rt_a_temp = x_solution' * mean(wk_return_d1_train, 2) - rf;
    risk_a = -1 * wk_return_d1_train' * x_solution - rho_solution;
    risk_a = max(0, risk_a);
    cvar_a_temp = rho_solution + 1/((i-1)*(1-theta)) * sum(risk_a);
    cvar_a(i-N/2) = cvar_a_temp;
    ratio_a_my = My_wk_rt_a_temp / cvar_a_temp;
    csr_a(i-N/2) = ratio_a_my;

    %equally weighted csr_p
    ew_wk_rt_temp = x_ew' * wk_return_d1_test - rf;
    ew_wk_rt(i-N/2) = ew_wk_rt_temp;   
    [~, cvar0_ew] = CVaR_p_cal(i-N/2, wk_return_d1, x_ew_all, theta);
    ratio_ew = sum(ew_wk_rt) / (i-N/2) / cvar0_ew;
    csr_ew(i-N/2) = ratio_ew;
    
    ew_ac_wk_rt_temp = ew_ac_wk_rt_temp + x_ew' * wk_return_d1_test;
    ew_ac_wk_rt(i-N/2) = ew_ac_wk_rt_temp;
    
    %index csr_p
    index_wk_rt_temp = wk_return_b1(i) - rf;
    index_wk_rt(i-N/2) = index_wk_rt_temp;

    f = [1 / ((i-1) * (1-theta)) * ones(1, i-1), 1];
    A = [-1*eye(i-1), -1*ones(i-1, 1);-1*eye(i-1), zeros(i-1, 1)];
    b = [wk_return_b1(1:i-1)'; zeros(i-1, 1)];
    solution_index = linprog(f, A, b);
    rho_index_solution = solution_index(i);

    risk_index_temp = -1 * wk_return_b1(:, i) - rho_solution;
    risk_index_temp = max(risk_index_temp, 0);
    risk_index(:, i-N/2) = risk_index_temp;
    cvar_index_temp = rho_index_solution + 1/((i-N/2)*(1-theta)) * sum(risk_index);
    cvar_index(i-N/2) = cvar_index_temp;
    ratio_index = sum(index_wk_rt) / (i-N/2) / cvar_index_temp;
    csr_index(i-N/2) = ratio_index;
    
    

end

ratio_p_all_yearly = ratio_p_all * (N/5.5)^0.5;
csr_a_yearly = csr_a * (N/5.5)^0.5;
csr_ew_yearly = csr_ew * (N/5.5)^0.5;
csr_index_yearly = csr_index * (N/5.5)^0.5;


%各个csr的最终值比较
ratio_p_all_yearly(end) 
csr_ew_yearly(end)       % 不太一样
csr_index_yearly(end)    % 还需要修正