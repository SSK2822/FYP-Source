function [c,ceq] = nonlinear_constrain(u, v0)
c =u'*v0*u-1;
ceq = [];