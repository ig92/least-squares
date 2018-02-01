% Conjugate Gradient Method
% Input:
%   - Q positive definite symmetric matrix
%   - q output matrix / vector
%   - x starting point
%   - eps precision
% Output:
%   - x minimising vector
function [x] = conjugate_gradient(Q, q, x, eps)
    g = Q * x - q;
    ng = norm(g);
    % the first direction is just the opposite of the gradient
    d = - g;
    while (ng > eps)
        ngs = ng^2;
        Qd = Q * d;
        alpha = (ngs / (d' * Qd));      % exact line search
        g = g + alpha * Qd;             % update of the gradient
        x = x + alpha * d;
        ng = norm(g);
        beta = ng^2 / ngs;              % Fletcher Reeves formula
        d = - g + (beta * d);           % update of the direction
    end
end