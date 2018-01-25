dims = [100 200 300 400 500 600 700 800 900 1000];

% testing linear complexity
A = [];
for j = 1 : size(dims,2)
    a = 0;
    for i = 1 : 10
        U = randn(dims(1,j), 100);
        tic; qr_mod(U); a = a + toc;
    end
    A = [A (a / 10)];
end
plot(A)

% testing quadratic complexity
A = [];
for j = 1 : size(dims,2)
    a = 0;
    for i = 1 : 10
        U = randn(1000, dims(1,j));
        tic; qr_mod(U); a = a + toc;
    end
    A = [A (a / 10)];
end
plot(A)