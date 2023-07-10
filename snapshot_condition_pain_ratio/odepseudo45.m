function obj = odepseudo45(y, miu)
obj = -(miu'*y); % 注意公式45原式是max，所以要加负号求最小值
end