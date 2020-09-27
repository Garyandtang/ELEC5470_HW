clear;clc;
n = 10
p = 8
X = zeros(n,p);
cMatrix = eye(n);
a = 2;
for i = 1:p
    noiseVec = randn(n,1);
    noiseVec = cMatrix*noiseVec;
    X(:,i) = noiseVec;
end
num = 100;
S = 1/n * (X')*X;
theta_1_norm = zeros(num,1);
for i = 1 : num
    cvx_begin
        variable theta(p,p) semidefinite
        minimize(trace(S*theta) - log_det(theta*S) + (i-1)*ones(1,p)*theta*ones(p,1))
    cvx_end
    theta_1_norm(i) = ones(1,p)*theta*ones(p,1);
end

a = 0 : num-1;
plot(a,theta_1_norm);
xlabel('${\alpha}$','Interpreter','latex','fontsize', 16)
ylabel('Quantity $||\Theta^*(\alpha)||$','Interpreter','latex','fontsize', 16)

