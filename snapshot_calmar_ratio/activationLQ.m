%定义激活函数 均匀分布的区间为[-1,0] y是高维的，所以g也是高维
function g = activationLQ(u)
    g = zeros(length(u),1);
    a = -1;
    b = 0;
    for i = 1:length(u)
        if u(i) < 0
            g(i) = -1;
        elseif u(i) == 0
            g(i) = a + (b-a)*rand(1);
        elseif u(i) > 0
            g(i) = 0;
        end
    end
end
