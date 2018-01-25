dim = 100;
U = randn(dim, dim);
U = U'*U;
u = randn(dim,1);
x0 = U * randn(dim,1);
w = conjugate_gradient(U,u,x0,10^-10)