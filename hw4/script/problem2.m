clc;clear;
rng('default')
rng(1);
n = 5;
Y = complex(normrnd(0, 1, n, n*2), normrnd(0, 1, n, n*2));
M = Y*Y';

% get maximum eigenvalue of M
[U, Sigma, V] = svd(M);
max_lambda = max(max(Sigma));

% MM algorithm
k = 1;
x = ones(n,1);
obj_value = [x'*M*x];
while 1
    x_ = x;
    y = -1*(max_lambda*eye(n) - M)*x;
    for i = 1 : n
        x(i) = -exp(angle(y(i))*1i);
    end
    
    obj_value = [obj_value; x'*M*x];
    k = k +1;
    if abs(x_ - x) <= 0.000001
        break
    end
end

plot([1:k],obj_value);
xlabel({'The numbber of iterations'},'Interpreter','latex','fontsize', 16);
ylabel({'Objective value'},'Interpreter','latex','fontsize', 16);
