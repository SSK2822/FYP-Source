clc,clear
format long;

load('NASDAQ100.mat'); %stock price matrix(939 weeks)
wk_return_d1 = Assets_Returns(1:596, :)';
wk_return_b1 = Index_Returns(1:596, :)';
[M, N] = size(wk_return_d1);
file1 = './NEW_NASDAQ100_Calmar/';
%file2 = './NEW_NASDAQ100_ESR_rho/';

rf_00_17 = load('rf_04_16a.txt',','); %risk free
rt_test = rf_00_17(N/2+1:N);

%Variable initial
theta = 0.95;
xt = zeros(M, 1);
xt_all = zeros(M, N/2);
x_ew = 1 / M * ones(M, 1);
ratio_a_all = zeros(1, N/2);
ratio_p_all = zeros(1, N/2);
nolinear_all = zeros(1, N/2);
rho_all = zeros(1, N/2);
My_wk_rt = zeros(1,N/2);

for i = (N/2+1):N %training output:xt，sor_p：只计算i=470这一周
    i
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %online setting 1-469ex-ante训练数据
    wk_return_d1_test = wk_return_d1(:, i); %只去second period的第一周
    rf = rf_00_17(i); %risk free
    miu = mean(wk_return_d1_train, 2); %miu
    if (miu < rf) %if no stock price greater than rf
        fprintf('do not trade in week %d\n', i)
        xt = zeros(M, 1); 
    else 
        [xt, ksi] = NNN_IG_linear(i, wk_return_d1, rf_00_17, miu); 
        %[xt, rho, nolinear] = NNN_IG(i, wk_return_d1, rf_00_17, miu); 
    end     
    xt_all(:, i-N/2) = xt;
    %nolinear_all(i-N/2) = nolinear; 
    %rho_all(i-N/2) = rho;
    
    %sor_p
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    %evar0 = EVaR_p_cal(i-N/2, wk_return_d1, xt_all, theta);
    % evar0 = EVaR_p(i-N/2, wk_return_d1, xt_all, rho, theta);
    %ratio_my = sum(My_wk_rt) / (i-N/2) / evar0;
    max_drawdown0 = Calmar_Var_p(i-N/2, wk_return_d1, xt_all);
    ratio_my = sum(My_wk_rt) / (i-N/2) / max_drawdown0;
    ratio_p_all(i-N/2) = ratio_my;
    
%{    
    %sor_a    
    if (i<N)
        rhot = rho;
        My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
        ratio_a_my = My_wk_rt_a_temp / EVaR_a(i, wk_return_d1, xt, rhot, theta);
        ratio_a_all(i-N/2) = ratio_a_my;
    end

    if (i == N)
        My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
        ratio_a_my = My_wk_rt_a_temp / EVaR_a_cal(i, wk_return_d1, xt, theta);
        ratio_a_all(i-N/2) = ratio_a_my;
    end
 %}   
    %保存每期的结果xt文件
    fid0 = fopen([file1, 'xt_', num2str(i), '.txt'], 'w'); %写入文件的命名
    for ii = 1:M
        fprintf(fid0, '%.6f\t', xt(ii, 1));
        fprintf(fid0, '\r\n');
    end
    fclose(fid0);
    
end

ratio_p_all* (N/18)^0.5; %update数据集
ratio_p_all* (N/5.5)^0.5; %or 数据集
rr = ratio_a_all* (N/18)^0.5;
rr(1:5)

ratiot_yearly = ratio_p_all*(N/26.25)^0.5;  %bsct DJIA数据集
ratiot_yearly = ratio_p_all* (N/11.5)^0.5;  %bsct NASDAQ100数据集
ratiot_yearly = ratio_p_all* (N/13.83)^0.5; %bsct FTSE100数据集
ratiot_yearly = ratio_p_all* (N/11.5)^0.5;