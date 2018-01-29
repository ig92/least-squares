blockdim = 100;
times = [];
n = 5;
dims = [];
while blockdim <= 1000
    time = 0;
    for i = 1 : n
        fprintf("solving qr with blockdim " + blockdim + " iteration " + i);
        tic; qr_mod_bw(A, blockdim);
        t = toc;
        time = time + toc;
        fprintf(" with total time " + t + "\n");
    end
    times = [times (time / n)];
    dims = [dims, blockdim];
    blockdim = blockdim + 100;
end
plot(dims,times);