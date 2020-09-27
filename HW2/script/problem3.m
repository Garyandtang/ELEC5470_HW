clear;clc;
n = randi([1,100],1);
m = randi([1,100],1);
X = randi([1,100],m,n);
y = randi([1,100],m,1);
cvx_begin
    variable x(n)
    minimize( (norm( y-X*x, 2 )) )
    subject to
        for i = 1:n-1
            x(i) <= x(i+1)
        end
cvx_end
solution = cvx_optval^2;
plot(x)
xlabel("index i",'fontsize', 16)
ylabel('Optimal solution element ${\beta_{i}}$','Interpreter','latex','fontsize', 16)
