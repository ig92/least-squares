function [w] = conjugate_gradient(Q, q, x0, eps)
    % Conjugate gradient method

    x = x0;
    
    % Gradient in the starting point x0
    g_new = grad(Q,q,x0);
    
    % At the first step, we take the direction exactly opposite
    % of the gradient in the starting point x0 (Gradient Descent)
    d_new = - g_new;
    
    while (norm(g_new) >= eps)
        % Exact line search (since strictly convex function)
        alpha = norm(g_new)^2 / (d_new' * Q * d_new);
        
        % Update the new point
        % Save the previous direction and gradient
        % used for the beta computation
        x = x + alpha * d_new;
        d_old = d_new;
        g_old = g_new;
        
        % Gradient in the new point x
        g_new = grad(Q,q,x);
        
        % Fletcher-Reevers formula for beta
        beta = norm(g_new)^2 / norm(g_old)^2;
        
        % The new direction
        d_new = - g_new + beta * d_old;
    end
    w = x;
end

function g = grad(Q, q, x)
    % Computes the gradient in point x of a quadratic function
    g = Q * x - q;
end