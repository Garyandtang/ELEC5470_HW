clear; clc;

% init parameter
randn('seed',1);
beta = zeros(10,1); beta(3) = 1; beta(5) = 7; beta(10) = 3;
% beta = ones(10,1);
n = 100; p = 10;
X = randn(n,p);
y = X*beta + 0.1*randn(n,1);
lambda = 0.2;
t = 20*ones(p,1);
x = [beta; t];

% solve with CVX
% probelm
% cvx_begin
%     variable xx(p)
%     minimize(sum_square_abs(y-X*xx)+lambda*norm(xx,1))
% cvx_end

% barrier mehtod
k = 1;
alpha = 5;
mu = 15;
tolerance = 10^(-6);

f = get_f(X, y, lambda, p, alpha, x);

f_list = [f];


dual_gap_all = [];
while (1)
    [dual_gap,x_, f_] = backtracking_newton(X,y,lambda,p,alpha,x,tolerance);
%     [x_, f_] = newton(X,y, lambda, p, alpha, x, tolerance);
    dual_gap_all = [dual_gap_all; dual_gap];
    x = x_;
    f = f_;
    f_list = [f_list; f];
    k = k + 1;
    if (2*p/alpha < tolerance)
        break
    end
    alpha = alpha*mu;
end

% draw
figure (1)
plot([1:k], f_list,'LineWidth',3)
xlabel("Barrier iteration i",'Interpreter','latex','fontsize', 40,'LineWidth',8)
ylabel(' Objective function value versus','Interpreter','latex','fontsize', 40,'LineWidth',8)
figure (2)
plot([1:size(dual_gap_all)], dual_gap_all,'LineWidth',3)
set(gca,'yscale','log')
xlabel("Newton iteration",'Interpreter','latex','fontsize', 40,'LineWidth',5)
ylabel('Duality gap bound','Interpreter','latex','fontsize', 40,'LineWidth',5)


%%% functions
function [x, f] = newton(X,y, lambda, p, alpha, x, tolerance)
    
    while (1)
        % computer newton step and decrement
        [g, H] = get_g_H(X,y, lambda, p, alpha, x);
    
        delta_x = -inv(H)*g;
        lambda_x_2 = g'*inv(H)*g;
        % stopping crterion: lambda(x)^2/2 <= tolerance
        if (lambda_x_2/2 <= tolerance)
            break;
        end
        % line search
        t = 1;
        b = 0.9; a = 0.1;
        newx = x + t* delta_x;
        while (get_f(X,y, lambda, p, alpha, newx) >= get_f(X,y, lambda, p, alpha, x) + a*t*g'*delta_x)
            t = b*t;
            newx = x + t* delta_x;
        end
        % update
        x = newx;
        f = get_f(X,y,lambda,p,alpha,x);
    end
    
end

function [dual_gap, newx,newf_value]=backtracking_newton(X,y,lambda,p,delta,x,error_tol)

    [g,H]=get_g_H(X,y,lambda,p,delta,x);
    deltax = -inv(H)*g;
    decrement_2 = g'*inv(H)*g;
    t = 0.01;
    newx = x;
    newf_value = get_f(X,y,lambda,p,delta,x);
    
    
   

    a = 0.1;
    b = 0.9;
    dual_gap = [];
    % 2.Stopping criterion. quit if λ2/2≤ error_tol.
    while (decrement_2/2>error_tol)
        
        % 3.Line search. Choose step size t by backtracking line search.
        newx = x+t*deltax;
        while (get_f(X,y,lambda,p,delta,newx) >= get_f(X,y,lambda,p,delta,x) +a*t*(g')*deltax)
            t = b*t;
            newx = x+t*deltax;
            
           
        end
        dual_gap = [dual_gap; 20/delta];
        % 4a.Update. x:=x+t∆xnt
        x = newx;
        newf_value = get_f(X,y,lambda,p,delta,x);
        

        
        % 1.Compute the Newton step and decrement.
        [g,H]=get_g_H(X,y,lambda,p,delta,x);
        deltax = -inv(H)*g;
        decrement_2 = g'*inv(H)*g;
    end
    
end

function [g, H] = get_g_H(X,y, lambda, p, alpha, x)
    beta = x(1:p);
    t = x(p+1:2*p);
    % computer the gradient
    g_b_0 = 2*X'*X*beta - 2*X'*y;
    g_b_1 = zeros(p,1);
    for i = 1: p
        g_b_1(i) = 1/alpha *(2*beta(i)/(t(i)^2-beta(i)^2));
    end
    g_b = g_b_0 + g_b_1;
    
    g_t_0 = lambda*ones(p,1);
    g_t_1 = zeros(p,1);
    for i = 1:10
        g_t_1 = -1/alpha*(2*t(i)/(t(i)^2-beta(i)^2));
    end
    g_t = g_t_0 +g_t_1;
    g = [g_b;g_t];
    
    % computer the hessian
    P1_v = zeros(p,1);
    for i=1:p
        P1_v(i) = 2*(t(i)^2+beta(i)^2)/(t(i)^2-beta(i)^2)^2/alpha;
    end
    P1 = diag(P1_v);

    P2_v = zeros(p,1);
    for i=1:p
        P2_v(i) = -4*(t(i)*beta(i))/(t(i)^2-beta(i)^2)^2/alpha;
    end
    P2 = diag(P2_v);

    H = [2*X'*X+P1, P2;P2,P1];
    
end

function f = get_f(X, y, lambda, p, alpha, x)
    beta = x(1:p);
    t = x(p+1:2*p);
    f = (y - X*beta)'*(y - X*beta) + lambda*ones(p,1)'*t;
    pha = 0;
    for i = 1:p
        pha = pha + log(t(i)-beta(i)) + log(beta(i)+t(i));
    end
    pha = -1/alpha*pha;
    f = f+ pha;
end

