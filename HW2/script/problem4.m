clear;clc;
% n = randi([1,100],1);
% p = randi([1,100],1);
n = 10;
p = 8;
X = zeros(n,p);
cMatrix = eye(n);

for i = 1:p
    noiseVec = randn(n,1);
    noiseVec = cMatrix*noiseVec;
    X(:,i) = noiseVec;
end
S = 1/n * (X')*X;

num = 100;
theta_1_norm = zeros(num,1);
theta_list = zeros(p,p, num);
for i = 1 : num
    cvx_begin
        variable theta(p,p) semidefinite
        minimize(trace(S*theta) - log_det(theta) + (i-1)*ones(1,p)*abs(theta)*ones(p,1))
    cvx_end
    theta_1_norm(i) = ones(1,p)*abs(theta)*ones(p,1);
    theta_list(:,:,i) = theta;
end

a = 0 : num-1;
plot(a,theta_1_norm);
xlabel('${\alpha}$','Interpreter','latex','fontsize', 16)
ylabel('Quantity $||\Theta^*(\alpha)||$','Interpreter','latex','fontsize', 16)

