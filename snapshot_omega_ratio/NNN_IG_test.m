function xk = NNN_IG(weeks, wk_return_d1, rf, miu) 
%weeks(>N/2),wk_return_d1: whole data，rb: whole data，miu: mean
wk_return_train = wk_return_d1(:, 1:weeks-1); %online
[M, ~] = size(wk_return_train); 
A = [ones(1,M),-1;miu',-rf(weeks)];
b = [0;1];
P1 = (A')*inv(A*A')*A;    % 直接计算新公式 此时P公式的A需要根据约束条件重写
P2 = 1 / M * ones(M, M); % 直接计算sharp ratio 的伪凸
I1 = eye(M+1);
I2 = eye(M); 
q1 = (A')*inv(A*A')*b;
q2 = 1 / M * ones(M, 1);



%variables initial
x0 = 1 / M * ones(M, 1); %initial as equally weighted
iter = 11; %10 times is enough
ga = zeros(iter, 1); %gamma_iteration

cov0 = x0' * cov(wk_return_train') * x0; %initially cov matrix
ga(1) = (miu' * x0 - rf(weeks)) / cov0; %initially gamma
gaa = ga(1);


iteration = 1;
while iteration < iter

    if (gaa < 0.000001) %condition: gamma>0
        n1 = length(miu(miu > rf(weeks)));
        for s = 1:M %更新x0
            if (miu(s) > rf(weeks))
                x0(s) = 1 / n1;
            end
        end 
        %if gamma<0, produce a new decision vector to make sure gamma>0
        cov0 = x0' * cov(wk_return_train') * x0;%renew cov0
        gaa = (miu' * x0 - rf(weeks)) / cov0; 
        ga(iteration) = gaa;  %renew gamma
    else 

        %ode hyperparameter 常微分方程参数
        step = 0.004;
        t_end = 0.08;
        tspan = 0:step:t_end; % 共21步
        if (iteration>4)
            step = 0.002;
            t_end = 0.01;
            tspan = 0:step:t_end; % 共6步
        end
        %ode solver(sometimes use ode15s if needed)
        [t,u] = ode45(@(t,u) odepseudoNN(u, P1, q1, I1, miu, wk_return_train, 0.0001), tspan, [x0;0.5]);
        [n, ~] = size(u);
        sr_plot = zeros(1, n);
        for i = 1:n
            u_plot = u(i, :);
            x_plot = max(0, u_plot);
            x_plot = x_plot';
            x_plot = x_plot(1:M)/x_plot(M); % y/ksi变换得到x
            My_wk_rt_a_temp = x_plot' * mean(wk_return_train, 2) - rf(weeks);
            ratio_a_my = My_wk_rt_a_temp / (x_plot'*cov(wk_return_train')*x_plot)^0.5; %as defined in sr_a
            sr_plot(i) = ratio_a_my;
        end
        t_plot = (iteration-1) * t_end * ones(n, 1) + t;
        if (iteration>4)
            t_plot = ((iteration-5) * t_end + 0.32) * ones(n, 1) + t;
        end
        figure(1);%k次迭代中sharp ratio的收敛情况
        sr_plot
        p1 = plot(t_plot, sr_plot*0.261249/0.251916, '-','linewidth', 1.5, 'markersize', 5,'color', 'r'); %迭代加权
        hold on
        u0 = u(end, :);
        x0 = max(0, u0);%relu
        x0 = x0';
        x0 = x0(1:M)/x0(M);
        cov0 = x0' * cov(wk_return_train') * x0;
        gaa = (miu' * x0 - rf(weeks)) / cov0; 
        ga(iteration) = gaa; 
    end
    %dsr_p(iteration) = DSR_p(weeks, wk_return_d1, x0, rf);
    iteration = iteration + 1;
end 
step_pseudo = 0.0002;
t_end_pseudo = 0.4;
tspan_pseudo = 0:step_pseudo:t_end_pseudo; 
[t_pseudo,u_pseudo] = ode45(@(t,u) odepseudo(u, P2, q2, I2, miu, wk_return_train, 0.0001, rf(weeks)), tspan_pseudo, 1 / M * ones(M, 1));
[m, ~] = size(u_pseudo);
sr_plot_pseudo = zeros(1, m);
for j = 1:m
    u_plot_pseudo = u_pseudo(j, :);
    x_plot_pseudo = max(0, u_plot_pseudo);
    x_plot_pseudo = x_plot_pseudo';
    My_wk_rt_a_temp = x_plot_pseudo' * mean(wk_return_train, 2) - rf(weeks);
    ratio_a_my = My_wk_rt_a_temp / (x_plot_pseudo'*cov(wk_return_train')*x_plot_pseudo)^0.5; %as defined in sr_a
    sr_plot_pseudo(j) = ratio_a_my;
end
t_plot_pseudo = t_pseudo;
linewidth_gca = 1.6;
fontsize_gca = 24;
fontsize_label = 28;
set(gca, 'linewidth', linewidth_gca, 'fontsize', fontsize_gca);
ylabel('\itSR', 'fontsize', fontsize_label)   ;
xlim([-0.002 0.38]);
ylim([0.15, 0.27])
title('SP 500'); 
hold on
p2 = plot(t_plot_pseudo, sr_plot_pseudo*0.261249/0.251916, '-.','linewidth', 1.5, 'markersize', 5,'color', 'b'); %直接算sharp ratio
set(gcf, 'unit', 'centimeters', 'position', [10, 10, 16, 10])
fontsize_legend = 20;
h = legend(p1, ['One step',sprintf('\n'), 'convex optimization']);
legend('boxoff')
ah=axes('position',get(gca,'position'),'visible','off');
h2 = legend(ah, p2, 'Pseudoconvex Optimization');
legend('boxoff')
set(h,'fontsize',fontsize_legend);
set(h2,'fontsize',fontsize_legend);
hold on

xk = x0;%output
end

