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
    D = [];
    R = [];
    i = 1;
    G = [];
    while (ng > eps)
        ngs = ng^2; % square of the norm of the gradient
        % this computation avoids computing twice the same thing
        Qd = Q * d; 
        alpha = ( ngs / (d' * Qd) );    % exact line search
        g = g + alpha * Qd;             % update of the gradient
        G = [G g];
        for j = 1 : i-1
            gamma = (g'*G(1:end,j))/(G(1:end,j)'*G(1:end,j));
            g = g - gamma * G(1:end,j)
        end
        i = i+1;
        ng = norm(g);
        beta = ng^2 / ngs;              % Fletcher Reeves formula
        D = [D (d ./ norm(d))];
        d = - g + (beta * d);           % update of the direction
        R = [R norm(D'*Q*d)];
    end
    w = x;
    plot(R)
end