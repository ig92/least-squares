dims = [100 200 300 400 500 600 700 800 900 1000];

% testing linear complexity
A = [];
B = [];
for j = 1 : size(dims,2)
    a = 0;
    b = 0;
    for i = 1 : 10
        fprintf("iteration " + i + " of dim " + dims(1,j) + "\n");
        % least squares with qr
        U = randn(dims(1,j), 100);
        u = randn(size(U,1), 1);
        tic; least_squares(U,u); a = a + toc;
        % least squares with cg
        u = U'*u;
        U = U' * U;
        x = U * randn(size(U,1),1);
        tic; conjugate_gradient(U, u, x, 10^-6); b = b + toc;
    end
    A = [A (a / 10)];
    B = [B (b / 10)];
end
plot(A,B)