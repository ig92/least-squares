dim = 100;
n = 5;
dims = [];
times = [];
while (dim < 1200)
    time = 0;
    for i = 1 : 10
        fprintf("solving dim " + dim + " iteration " + i);
        AA = randn(dim,dim);
        AA = AA' * AA;
        aa = randn(dim,1);
        x0 = AA * randn(dim,1);
        tic; conjugate_gradient(AA,aa,x0,10^-10); 
        t = toc;
        time = time + toc;
        fprintf(" in time " + t + "\n");
    end
    time = toc;
    times = [times time / n];
    dims = [dims dim];
    dim = dim + 50;
    
end
plot(dims, times)