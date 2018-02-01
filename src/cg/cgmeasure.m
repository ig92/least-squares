% cgmeasure is a simple batch that measures the time of
% conjugate gradient method
% Input:
%   - dimStart starting dimension of Q square matrix to test
%   - dimEnd the last dimension of Q square matrix to test
%   - dimIncrease step size of the increase of the dimension
%   - n number of iterations
%   - eps precision
% Output:
%   - dims all the dimensions that have been tested
%   - times the times of corresponding dimensions in dims
function [dims, times] = cgmeasure(dimStart, dimEnd, dimIncrease, n, eps)
    dim = dimStart;
    dims = [];
    times = [];
    while (dim <= dimEnd)
        time = 0;
        for i = 1 : n
            fprintf("solving dim " + dim + " iteration " + i);
            Q = randn(dim,dim);
            Q = Q' * Q;
            q = randn(dim,1);
            x0 = Q * randn(dim,1);
            tic; conjugate_gradient(Q,q,x0,eps); 
            t = toc;
            time = time + toc;
            fprintf(" in time " + t + "\n");
        end
        time = toc;
        times = [times time / n];
        dims = [dims dim];
        dim = dim + dimIncrease;
    end
end