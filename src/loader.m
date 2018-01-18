cd /code/matlab/least-squares/src/
X = csvread('../test/ml_extended.csv');
b = csvread('../test/ml-output.csv');
Q = X'*X;
b1 = b(1:end,1);
b2 = b(1:end,2);
q1 = X'*(b1);
q2 = X'*(b2);
x0 = randn(12,1);
eps = 10^-10;