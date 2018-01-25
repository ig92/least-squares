function [w] = conjugate_gradient(Q, q, x0, eps)
% conjugate gradient method
% Q is the symmetric square matrix positive definite
% q is a vector
% x0 is the starting point
% eps is the required precision

    x = x0;
    g = Q * x - q;
    ng = norm(g);
    % the first direction is the opposite of the gradient
    d = - g;
    while (ng > eps)
        ngs = ng^2; % square of the norm of the gradient
        % this computation avoids computing twice the same thing
        Qd = Q*d;
        alpha = ( ngs / (d' * Qd) );    % exact line search
        x = x + alpha * d;
        g = g + alpha * Qd;             % update of the gradient
        ng = norm(g);
        beta = ng^2 / ngs;              % Fletcher Reeves formula
        d = - g + (beta * d);           % update of the direction
    end
    w = x;
end