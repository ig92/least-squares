function [Q,R] = qr_mod_bw(A,blockdim)
    [m,n] = size(A);

    % compute upper triangular R
    H = zeros(m,n);

    blockstart = 1;
    blockend = blockdim;
    hh = 1;

    while blockstart < n
        % if the last block is smaller than blockdim
        if n < blockend
            blockend = n;
        end

        % take the new block
        B = A(1:end,blockstart:blockend);

        % apply householders to the new block starting from the second one
        if (blockstart > 1)
            B = apply_householders(H,B);
            % take upper part since it is already a peace of the final R
            A(1:hh-1,blockstart:blockend) = B(1:hh-1,1:end);
            % take lower part for further processing
            B = B(hh:end,1:end);
        end

        % create new householders for the new block
        for i = 1 : size(B,2)
            v = hh_reflector(B(i:end,i));
            H(hh:end,hh) = v;
            B(i:end,i:end) = B(i:end,i:end) - 2 * v * (v' * B(i:end,i:end));
            B(i+1:end,i) = 0;
            hh = hh + 1;
        end

        % copy the computer upper triangular on B
        A(hh-blockdim:end,blockstart:blockend) = 0;
        A(hh-blockdim:hh-1,blockstart:blockend) = B(1:blockdim,1:end);

        % update block delimiters
        blockstart = blockend + 1;
        blockend = blockend + blockdim;
    end

    % compute orthogonal Q from v_i's
    Q = eye(m,n);
    bs = (n - blockdim) + 1;
    be = n;
    while bs > 0
        if (bs < 1)
            bs = 1;
        end
        B = Q(1:end,bs:be);
        for i = be : -1 : 1
            v = H(i:end,i);
            B(i:end,1:end) = B(i:end,1:end) - 2 * (v * (v' * B(i:end,1:end)));
        end
        Q(1:end,bs:be) = B;
        be = bs - 1;
        bs = bs - blockdim;
    end

    % drop zeros
    R = A(1:n,1:end);
end

function B = apply_householders(H,B)
    m = size(B,1);
    for i = 1 : (size(H,2))
        v = H(i:(m+1)-1,i);
        B(i:end,1:end) = B(i:end,1:end) - 2 * v * (v' * B(i:end,1:end));
    end
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
end