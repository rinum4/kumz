fun = @(x) (x(1,1) - 0.2)^2 + ...
    (x(1,2) - 1.7)^2 + (x(2,1) - 5.1)^2 - x(2,2);
x = ga(fun,4,[],[],[],[],[],[],[], [2 3]);