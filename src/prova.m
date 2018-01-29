[Q,R] = qr(A,0)
[Q1,R1] = qr(A(1:end,1:2));
Q1'*A