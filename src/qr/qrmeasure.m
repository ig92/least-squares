% qrmeasure is a simple batch that measure the time of qr factorisation.
% It is divided in two parts, in order to see the linearity in the number of
% of rows and the quadratica complessity in the number of rows.
% Input:
%   - startDim is the starting dimension
%   - endDim is the last dimension we check
%   - dimIncrease is the increasing step
%   - n is the number of times we iterate each matrix
% Output:
%   - dimsLin is a vector that contains the dimensions used for linear case
%   - timesLin is a vector with the times employed for linear case
%   - dimsQuad is a vector that contains the dimensions used for quadratic
%   case
%   - timesQuad is a vector with the times employed for quadratic case
function [dimsLin,timesLin,dimsQuad,timesQuad] = qrmeasure(startDim, endDim, dimIncrease, n)
    % linear testing
    dim = startDim;
    timesLin = [];
    dimsLin = [];
    while dim <= endDim
        time = 0;
        for i = 1 : n
            fprintf("linear with dimension " + dim + " iteration " + i);
            U = randn(dim, 100);
            u = randn(size(U,1), 1);
            tic; qr_least_squares(U,u); time = time + toc;
            t = toc;
            fprintf(" solved in time " + t + "\n");
        end
        dimsLin = [dimsLin dim];
        timesLin = [timesLin (time / n)];
        dim = dim + dimIncrease;
    end

    % quadratic testing
    dim = startDim;
    timesQuad = [];
    dimsQuad = [];
    while dim <= endDim
        time = 0;
        for i = 1 : n
            fprintf("quadratic with dimension " + dim + " iteration " + i);
            U = randn(1000, dim);
            u = randn(size(U,1), 1);
            tic; qr_least_squares(U,u); time = time + toc;
            t = toc;
            fprintf(" solved in time " + t + "\n");
        end
        dimsQuad = [dimsQuad dim];
        timesQuad = [timesQuad (time / n)];
        dim = dim + dimIncrease;
    end
end