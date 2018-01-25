
a = 0;
for i = 1 : 10
    U = randn(50000,100);
    tic; qr_mod(U); a = a + toc;
end

a / 10