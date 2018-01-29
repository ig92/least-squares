function [m, w] = least_squares(X, b)
% least_squares solves the least squares minimisation problem
% returns m which is the minimum value of the function of the
% least squares problem and w which is the vector that minimises
% the problem with input X and b

[Q1,R1] = qr_mod(X);    % QR factorise X
c = Q1' * b;            % intermediate step
w = R1 \ c;             % solved by backsubstitution
m = norm(X*w - b);      % the solution of least square