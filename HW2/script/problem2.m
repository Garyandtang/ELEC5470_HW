sigma = [1 0.0015 -0.02;
         0.0015 1 -0.1;
         -0.02 -0.1 1];
mu = [0.001,0.05,0.005]';
n = 3;

lambda  = zeros(11,1);
for i = 1:11
    lambda(i) = (i-1)*10^-4;
end

w_optimal = zeros(11,3);
for i = 1:11
    cvx_begin
        variable w(n)
        minimize(w'*sigma*w - lambda(i)*mu'*w)
        subject to
            w >= 0
            w'*ones(n,1) == 1
    cvx_end 
    w_optimal(i,:) = w';
end
expected_return = w_optimal*mu;
volatility = zeros(11,1);
for i = 1:11
    volatility(i)=sqrt(w_optimal(i,:)*sigma*w_optimal(i,:)');
end
<<<<<<< HEAD
plot(expected_return,volatility)
=======
plot(expected_return,(volatility))
xlabel('Portfolio Expected Return', 'fontsize',14)
ylabel('Portfolio Volatility','fontsize', 14)
>>>>>>> a1ce2feacf7e10f011fb3fddcab1d9e573f742e8
