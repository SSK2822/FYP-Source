clc,clear
%{
different dataset use:
FTSE:
file = './FTSE_DSR55/';
d1 = load('wk_price_FTSE', ',')'; 
b1 = load('b1_FTSE.txt', ','); 
HDAX:
file = './HDAX_DSR55/';
d1 = load('wk_price_HDAX', ',')'; 
b1 = load('b1_HDAX.txt', ','); 
HSCI:
file = './HSCI_DSR55/';
d1 = load('wk_price_HSCI', ',')';
b1 = load('b1_HSCI.txt', ','); 
SP500:
file = './SP_DSR55/';
d1 = load('wk_price_SP500', ',')'; 
b1 = load('SP500_historical_index', ',')';  
ORHDAX:
d1 = load('ORalldata_DAX.txt', ','); %price矩阵(939天数据)
b1 = load('ORb1_DAX.txt', ','); %股指 
ORFTSE:
d1 = load('ORalldata_FTSE.txt', ','); %price矩阵(939天数据)
b1 = load('ORb1_FTSE.txt', ','); %股指 
ORHSCI:
d1 = load('ORalldata_HS.txt', ','); %price矩阵(939天数据)
b1 = load('ORb1_HS.txt', ','); %股指 
ORSP100:
d1 = load('ORalldata_SP100.txt', ','); %price矩阵(939天数据)
b1 = load('ORb1_SP100.txt', ','); %股指 
%}

clc,clear

file = './BCST_SP500_SR/';
load('SP500.mat'); %stock price matrix(939 weeks)
wk_return_d1 = Assets_Returns(1:594, :)';
wk_return_b1 = Index_Returns(1:594, :)';
rf_00_17=load('rf_04_16b.txt',','); %risk free


%{
d1 = load('wk_price_FTSE', ',')'; 
b1 = load('b1_FTSE.txt', ','); 
file = './UPDATE_FTSE_SR/';

% rf_00_17=load('rf_92_97.txt',','); %risk-free
rf_00_17 = load('rf_00_17.txt',','); %risk-free

[wk_return_d1, ~] = price2ret(d1', [], 'Periodic'); %simple return
wk_return_d1 = wk_return_d1';
[wk_return_b1, ~] = price2ret(b1, [], 'Periodic'); %simple return

%}

[M, N] = size(wk_return_d1);
%Variable initial
xt = zeros(M, 1);
xt_all = zeros(M, N/2);
x_ew = 1 / M * ones(M, 1);
ratio_a_all = zeros(1, N/2);
ratio_p_all = zeros(1, N/2);
nolinear_all = zeros(1, N/2);
My_wk_rt = zeros(1,N/2);


for i = (N/2+1):N %training output:xt，sor_p：
    i
    wk_return_d1_train = wk_return_d1(:, 1:i-1); %online setting 1-469ex-ante训练数据
    wk_return_d1_test = wk_return_d1(:, i); %只去second period的第一周
    rf = rf_00_17(i); %risk free
    miu = mean(wk_return_d1_train, 2); %miu
    if (miu < rf) %if no stock price greater than rf
        fprintf('do not trade in week %d\n', i)
        xt = zeros(M, 1); 
    else 
        [xt, ratio, nolinear] = NNN_IG(i, wk_return_d1, rf_00_17, miu); 
    end
    xt_all(:, i-N/2) = xt;
    nolinear_all(i-N/2) = nolinear; 
    ratio_a_all(i-N/2) = ratio; % 先验
    
    %sr_p and cumulative return 后验
    My_wk_rt_temp = xt' * wk_return_d1_test - rf; %used to calculate denominater
    My_wk_rt(i-N/2) = My_wk_rt_temp;
    V_my = cov(My_wk_rt(1:i-N/2));%as defined in sr_p
    ratio_my = sum(My_wk_rt) / (i-N/2) / (V_my ^ 0.5);%本期sr值
    ratio_p_all(i-N/2) = ratio_my;%公式33的后验
    
    %保存每期的结果xt文件
    fid0 = fopen([file, 'xt_', num2str(i), '.txt'], 'w'); %写入文件的命名
    for ii = 1:M   
        fprintf(fid0, '%.6f\t', xt(ii, 1));
        fprintf(fid0, '\r\n');
    end
    fclose(fid0);
end

ratiot_yearly = ratio_p_all* (N/18)^0.5;       %update数据集
ratiot_yearly = ratio_p_all* (N/5.5)^0.5;    %or 数据集
ratiot_yearly = ratio_p_all*(N/26.25)^0.5;  %bsct DJIA数据集
ratiot_yearly = ratio_p_all* (N/11.5)^0.5;  %bsct NASDAC100数据集
ratiot_yearly = ratio_p_all* (N/13.83)^0.5; %bsct FTSE100数据集
ratiot_yearly = ratio_p_all* (N/11.5)^0.5;



