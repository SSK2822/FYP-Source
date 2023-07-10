function xk = NN(weeks, wk_return_d1, rf, miu) 
%weeks(>N/2),wk_return_d1: whole data，rb: whole data，miu: mean
wk_return_train = wk_return_d1(:, 1:weeks-1); %online
[M, ~] = size(wk_return_train); 
P = 1 / M * ones(M, M);
I = eye(M);
q = 1 / M * ones(M, 1);


%variables initial
x0 = 1 / M * ones(M, 1); %initial as equally weighted
iter = 10; %10 times is enough
ga = zeros(iter, 1); %gamma_iteration

cov0 = x0' * cov(wk_return_train') * x0; %initially cov matrix
ga(1) = (miu' * x0 - rf(weeks)) / cov0; %initially gamma
gaa = ga(1);

step = 0.0002;
t_end = 0.2;
tspan = 0:step:t_end; 
[t,u] = ode45(@(t,u) odepseudo(u, P, q, I, miu, gaa, wk_return_train, 0.0001, rf(weeks)), tspan, x0);
[n, ~] = size(u);
sr_plot = zeros(1, n);
for i = 1:n
    u_plot = u(i, :);
    x_plot = max(0, u_plot);
    x_plot = x_plot';
    My_wk_rt_a_temp = x_plot' * mean(wk_return_train, 2) - rf(weeks);
    ratio_a_my = My_wk_rt_a_temp / (x_plot'*cov(wk_return_train')*x_plot)^0.5; %as defined in sr_a
    sr_plot(i) = ratio_a_my;
end
t_plot = t;
figure(1);%k次迭代中gamma的收敛情况
plot(t_plot, sr_plot, '-.','linewidth', 1.5, 'markersize', 5,'color', 'b');
hold on
u0 = u(end, :);
x0 = max(0, u0);%relu
x0 = x0';
xk = x0;%output
end