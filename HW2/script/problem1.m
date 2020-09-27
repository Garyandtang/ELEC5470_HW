% problem 1
sigma = [1 0.0015 -0.02;
         0.0015 1 -0.1;
         -0.02 -0.1 1];
b = [ 0.1594, 0.0126, 0.8282]';
n = 3;

cvx_begin
    variable x(n)
    minimize(0.5 * x' * sigma * x - b'*log(x))
cvx_end

w = x/(ones(1,n)*x)

