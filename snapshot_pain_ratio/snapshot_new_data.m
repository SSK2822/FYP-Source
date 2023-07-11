clc,clear
format long;
%{
file = './BCST_DJIA_SOR/';
load('DowJones.mat'); %stock price matrix(939 weeks)
wk_return_d1 = Assets_Returns(2:1363, :)';
wk_return_b1 = Index_Returns(2:1363, :)';

rf_00_17=load('rf_90_16.txt',','); %risk free DJIA
%rf_00_17=load('rf_02_16.txt',','); %risk free FTSE
%}

file = './BCST_NASDAQ100_SOR/';
load('NASDAQ100.mat'); %stock price matrix(939 weeks)
wk_return_d1 = Assets_Returns(1:596, :)';
wk_return_b1 = Index_Returns(1:596, :)';
rf_00_17=load('rf_04_16a.txt',','); %risk free

%rf_00_17 = load('rf_00_17.txt',','); %risk-free update数据集
%rf_00_17 = load('rf_92_97.txt',','); %risk-free or数据集

[M, N] = size(wk_return_d1);

%Variable initial
xt = zeros(M, 1);
xt_all = zeros(M, N/2);
x_ew = 1 / M * ones(M, 1);
ratio_a_all = zeros(1, N/2);
ratio_p_all = zeros(1, N/2);
nolinear_all = zeros(1, N/2);
My_wk_rt = zeros(1,N/2);

%for i = (N/2+1):N %training output:xt，sor_p：只计算i=470这一周
for i = (N/2+1):(N/2+5) %training output:xt，sor_p：只计算i=470这一周
    i
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %online setting 1-469ex-ante训练数据
    wk_return_d1_test = wk_return_d1(:, i); %只去second period的第一周
    rf = rf_00_17(i); %risk free
    miu = mean(wk_return_d1_train, 2); %miu
    if (miu < rf) %if no stock price greater than rf
        fprintf('do not trade in week %d\n', i)
        xt = zeros(M, 1); 
    else 
        [xt, nolinear] = NNN_IG(i, wk_return_d1, rf_00_17, miu); 
    end     
    xt_all(:, i-N/2) = xt;
    nolinear_all(i-N/2) = nolinear; 
    
    %sor_p
    My_wk_rt_temp = xt' * wk_return_d1_test - rf;
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    % dsr0 = DSR_p(i, wk_return_d1, xt_all, rf_00_17);
    dsr0 = DSR_p_or(i-N/2, wk_return_d1, xt, rf_00_17);
    ratio_my = sum(My_wk_rt) / (i-N/2) / (dsr0 ^ 0.5);
    ratio_p_all(i-N/2) = ratio_my;
    
    %sor_a
    My_wk_rt_a_temp = xt' * mean(wk_return_d1_train, 2) - rf;
    ratio_a_my = My_wk_rt_a_temp / (DSR_a(i, wk_return_d1, xt, rf_00_17))^0.5;
    ratio_a_all(i-N/2) = ratio_a_my;
    
    %保存每期的结果xt文件
    fid0 = fopen([file, 'xt_', num2str(i), '.txt'], 'w'); %写入文件的命名
    for ii = 1:M
        fprintf(fid0, '%.10f\t', xt(ii, 1));
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