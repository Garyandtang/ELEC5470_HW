clear;clc;
n = randi([1,100],1);
p = randi([1,n-1],1);
X = zeros(n,p);
cMatrix = eye(n);
a = 0;
for i = 1:p
    noiseVec = randn(n,1);
    noiseVec = cMatrix*noiseVec;
    X(:,i) = noiseVec;
end

S = 1/n * (X')*X;

cvx_begin
    variable theta(p,p) semidefinite
    minimize(trace(S*theta) - log_det(theta*S)+ a*ones(1,p)*theta*ones(p,1))
    subject to
        theta <In> semidefinite(p);
cvx_end
        