%定义常微分方程
function dudt = odepseudoGSM(u, P1, q1, I, miu, wk_return_d1_train, ep1)%ep是epsilon

%syms x;
%PNN  (应该用P1 q1 I1矩阵 因为约束条件变了)
x = max(0,u);%relu（u）
g = (I - P1) * x + q1;%注意是57*1 M+1维

V = cov(wk_return_d1_train'); % 56*56
V = [[V;zeros(1,length(V))],zeros(length(V)+1,1)]; % 多了一个变量，维数+1，要在后面补上0
dfx = V * g /(g'*V*g)^0.5;  
dudt = 1 / ep * [-P1 * x - (I - P1) * (u - x + dfx) + q1]; 

%one-layer recurrent NN 有限时间收敛
%activation function
%y = u(1:M);
%ita = u(M+1);
%equ1 = ones(1,M)*y - ita;
%equ2 = miu'*y - rf*ita -1;
%V = cov(wk_return_d1_train');
%dydt = 1 / ep1 *[-V * y /(y'*V*y)^0.5 - alpha*ones(M,1)*activationEQU(equ1) - alpha*miu*activationEQU(equ2) - beta*activationLQ(y)];
%ditadt = 1 / ep2 *[alpha*activationEQU(equ1) + alpha*rf(weeks)*activationEQU(equ2) - beta*activationLQ(ita)];
 

end