function f = make_interpolant(x, y)
x = x(:); 
y = y(:);

[xs, idx] = sort(x);
ys       = y(idx);
[xu, ia] = unique(xs, 'stable');
yu       = ys(ia);

f = @(xq) interp1(xu, yu, xq, 'pchip');
end