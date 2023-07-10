clear;
clc;
close all;
% tic
tStart = cputime;
n=500;
x = sdpvar(n,1);
%y = binvar(n,1);
Con = [
sum(max(0,x).^2)<=1;
];
ops = sdpsettings('verbose',2,'solver','cplex');%,'maxtime',5555
ops.bnb.maxtime=3600*6;
z = sum((x+1).^2);
result = optimize(Con,z,ops);%
if result.problem ==0
x = value(x)
value(z)
else
disp('Something went wrong');
end
tEnd = cputime - tStart;