%定义激活函数 均匀分布的区间为[-1,1]
function g = activationEQU(u)
    a = -1;
    b = 1;
    if u < 0
        g = -1;
    elseif u == 0
        g = a + (b-a)*rand(1);
    elseif u > 0
        g = 1;
    end
end
