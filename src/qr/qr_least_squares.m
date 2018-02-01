% Least Squares with QR Factorisation
% Input:
%   - X a matrix
%   - b a vector
% Output:
%   - m minimum value
%   - w minimising vector / matrix
function [m,w] = qr_least_squares(X, b)
    [Q1,R1] = qr_mod(X);    % QR factorise X
    c = Q1' * b;            % intermediate step
    w = R1 \ c;             % solved by backsubstitution
    m = norm(X*w - b);      % the solution of least square
end