%定义常微分方程 直接计算sharp ratio 伪凸
function dudt = odepseudo(u, P, q, I, miu, wk_return_d1_train, ep, rf)%ep是epsilon

%syms x;

x = max(0,u);%relu（u）

g = (I - P) * x + q;%
%fprintf('good g %f\n',g);%465个,48个的和为1，全为正数
%dfx=(2*V*g*(g'*miu)-g'*V*g*miu)'/(g'*miu).^2; %df(x) change ~~~

V = cov(wk_return_d1_train');
dfx = -1 * (miu*(g'*V*g)^0.5-(miu'*g-rf)*V*g/(g'*V*g)^0.5) / (g'*V*g);


dudt = 1 / ep * [-P * x - (I - P) * (u - x + dfx) + q];  %equation (16)



end