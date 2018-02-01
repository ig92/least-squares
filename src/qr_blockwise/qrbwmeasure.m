% qrbwmeasure is a simple batch that evaluates
% the improvements obtanied with the block wise version of qr.
% At each iteration we execute the qr block wise with a random
% matrix and take the time employed.
% Input:
%   - r is the number of rows
%   - c is the number of columns
%   - n is the number of times we execute block wise qr with the same block
%   dimension
%   - blockStart is the block dimension employed at the beginning
%   - blockEnd is the bigger dimension we test
%   - blockIncrement is the increment of the block dim
% Output (used to create the graph):
%   - dims is a vector with the dimensions used 
%   - times is a vector with the times computed
function [dims,times] = qrbwmeasure(r,c,n,blockStart,blockEnd,blockIncrement)
    blockdim = blockStart;
    times = [];
    dims = [];
    while blockdim <= blockEnd
        time = 0;
        for i = 1 : n
            fprintf("solving qr with blockdim " + blockdim + " iteration " + i);
            A = randn(r,c);
            tic; qr_mod_bw(A,blockdim);
            t = toc;
            time = time + toc;
            fprintf(" with total time " + t + "\n");
        end
        times = [times (time / n)];
        dims = [dims, blockdim];
        blockdim = blockdim + blockIncrement;
    end
end