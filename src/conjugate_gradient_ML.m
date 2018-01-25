clear
function [m,w] = conjugate_gradient_ML(X,b)
    Q = X'*X;
    b1 = b(1:end,1);
    b2 = b(1:end,2);
    q1 = X'*(b1);
    q2 = X'*(b2);
    x0 = Q * randn(11,1);
    eps = 10^-10;
    w1 = conjugate_gradient(Q, q1, x0, eps)
    w2 = conjugate_gradient(Q, q2, x0, eps)
    w = [w1 w2];
    m = norm(X*w - b);
end