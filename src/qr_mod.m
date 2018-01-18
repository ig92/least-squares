function [Q, R] = qr_mod(A)
% computes QR of A

[m,n] = size(A);
min_ = min(m,n); % used for handling both short fat and tall thin cases

% compute upper triangular R
R = [A; zeros(1,n)];
for i = 1 : min_
    [v,R(i,i)] = hh_reflector(R(i:end-1,i));
    R(i+1:end,i) = v;
    R(i:end-1,i+1:end) = R(i:end-1,i+1:end) - 2 * v * (v' * R(i:end-1,i+1:end));
end

% compute orthogonal Q from v_i's
Q = eye(m,min_);
for i = min_ : -1 : 1
    v = R(i+1:end,i);
    Q(i:m,i:min_) = Q(i:m,i:min_) - 2 * (v * (v' * Q(i:m,i:min_)));
end

% remove v_i's from R
R = R(1:min_,1:end);
for i = 1 : min_ - 1
    R(i+1:end,i) = 0;
end

function [v,n] = hh_reflector(x)
% computes x - norm(x) * eye(size(x,1), 1)
if (x(1) ~= 0) % handles the case of the first entry 0
    n = - sign(x(1)) * norm(x);
else
    n = - norm(x);
end
v = x;
v(1) = v(1) - n;
if (norm(v) ~= 0) % handles the case of the zero vector
    v = v / norm(v);
end