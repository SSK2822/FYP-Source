clc,clear
format long;

d1 = load('wk_price_HSCI', ',')'; %price矩阵(939天数据)
b1 = load('b1_HSCI.txt', ','); %股指 
file1 = './UPDATE_HSCI_Martin/';

rf_00_17 = load('rf_00_17.txt',','); %risk-free or数据集

[wk_return_d1, ~] = price2ret(d1', [], 'Periodic'); %simple return
wk_return_d1 = wk_return_d1';
[wk_return_b1, ~] = price2ret(b1, [], 'Periodic'); %simple return
[M, N] = size(wk_return_d1);

%Variable initial
theta = 0.95;
xt = zeros(M, 1);
xt_all = zeros(M, N/2);
x_ew = 1 / M * ones(M, 1);
ratio_a_all = zeros(1, N/2);
ratio_p_all = zeros(1, N/2);
nolinear_all = zeros(1, N/2);
rho_all = zeros(1, N/2);
ksi_all = zeros(M, N/2);
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
        [xt, ksi, nolinear] = NNN_IG(i, wk_return_d1, rf_00_17, miu); 
    end     
    xt_all(:, i-N/2) = xt;
    % nolinear_all(i-N/2) = nolinear; 
    % rho_all(i-N/2) = rho;
    % ksi_all(:, i-N/2) = ksi;
    
    %martin_ratio_p
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    
    martin_drawdown0 = Martin_Var_p(i-N/2, wk_return_d1, xt_all);
    ratio_my = sum(My_wk_rt) / (i-N/2) / (martin_drawdown0)^0.5;
    ratio_p_all(i-N/2) = ratio_my;
   
    %martin_ratio_a
    My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
    ratio_a_my = My_wk_rt_a_temp / (Martin_Var_a(i, wk_return_d1, xt))^0.5;
    ratio_a_all(i-N/2) = ratio_a_my;

    
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