% Conjugate Gradient for ML Cup
% Input:
%   - X the \hat{X} matrix augmented from ML Cup dataset
%   - b the output vector: y_1 or y_2
%   - eps precision
% Output:
%   - m minimum value
%   - w minimising vector / matrix
function [m,w] = conjugate_gradient_ML(X,b,eps)
    Q = X'*X;
    b1 = b(1:end,1);
    b2 = b(1:end,2);
    q1 = X'*(b1);
    q2 = X'*(b2);
    x0 = Q * randn(size(Q,1),1);
    w1 = conjugate_gradient(Q, q1, x0, eps);
    w2 = conjugate_gradient(Q, q2, x0, eps);
    w = [w1 w2];
    m = norm(X*w - b);
end