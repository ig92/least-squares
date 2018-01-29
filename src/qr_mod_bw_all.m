function [Q, R] = qr_mod_bw(A, blockdim)
% computes QR of A

[m,n] = size(A);

% compute upper triangular R
R = A;
H = zeros(m,n);

blockstart = 1;
blockend = blockdim;

while blockstart < n
    if n < blockend
        blockend = n;
    end
    if (blockstart > 1)
        for j = 1 : blockstart-1
            v = H(j:end,j);
            R(j:end,blockstart:blockend) = R(j:end,blockstart:blockend) - 2 * v * (v' * R(j:end,blockstart:blockend));
        end
    end
    for i = blockstart : blockend
        [v,R(i,i)] = hh_reflector(R(i:end,i));
        H(1+(m-size(v,1)):end,i) = v;
        R(i+1:end,i) = zeros(size(v,1)-1, size(v,2));
        R(i:end,i+1:blockend) = R(i:end,i+1:blockend) - 2 * v * (v' * R(i:end,i+1:blockend));
    end
    blockstart = blockend + 1;
    blockend = blockend + blockdim;
end

% compute orthogonal Q from v_i's
Q = eye(m,n);
bs = (n - blockdim) + 1;
be = n;
j = n;
while bs > 0
    if (bs < 1)
        bs = 1;
    end
    for i = be : -1 : bs
        v = H(i:end,i);
        Q(j:m,j:n) = Q(j:m,j:n) - 2 * (v * (v' * Q(j:m,j:n)));
        j = j-1;
    end
    be = bs - 1;
    bs = bs - blockdim;
end

% remove v_i's from R
R = R(1:n,1:end);

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