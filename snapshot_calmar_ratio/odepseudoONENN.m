%定义常微分方程
function dy_dita_dt = odepseudoONENN(u, miu, wk_return_d1_train, rf, weeks, eps1, eps2, alpha, beta)%ep是epsilon

[M,~] = size(wk_return_d1_train);
dy_dita_dt = zeros(M+1, 1);

%syms x;
%one-layer recurrent NN 有限时间收敛
%activation function
y = u(1:M);
ita = u(M+1);
equ1 = ones(1,M)*y - ita;
equ2 = miu'*y - rf(weeks)*ita -1;
V = cov(wk_return_d1_train');
dy_dita_dt(1:M) = 1 / eps1 * (-(V * y )/(y'*V*y)^0.5 - alpha*ones(M,1)*activationEQU(equ1) - alpha*miu*activationEQU(equ2) - beta*activationLQ(y));
dy_dita_dt(M+1) = 1 / eps2 * (alpha*activationEQU(equ1) + alpha*rf(weeks)*activationEQU(equ2) - beta*activationLQ(ita));

 
end