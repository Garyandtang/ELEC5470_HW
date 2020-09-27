## ELEC5470 HW2

TANG Jiawei 20672550

### Problem1

$$\underset{\boldsymbol{x}>0}{\operatorname{minimize}} \quad \frac{1}{2} \boldsymbol{x}^{\top} \Sigma \boldsymbol{x}-\boldsymbol{b}^{\top} \log (\boldsymbol{x})$$

**Analysis**

As both quadratic function,$\frac{1}{2} \boldsymbol{x}^{\top} \Sigma \boldsymbol{x}$ , and minus log function, $-\log (\boldsymbol{x})$, are convex function in $\mathbb{R}^n_+$. This RPP problem is convex.

**Code**

```matlab
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
```

with answer: $x^* =[0.4085, 0.1671, 0.9226]^T$, $w = [0.2726, 0.1115, 0.6158]^T$

如何 subject to x >0 ????????????????????????????????????????????????

### Problem2

$$\begin{array}{ll}
\underset{w}{\operatorname{minimize}} & \boldsymbol{w}^{\top} \boldsymbol{\Sigma} \boldsymbol{w}-\lambda \boldsymbol{\mu}^{\top} \boldsymbol{w} \\
\text { subject to } & \boldsymbol{w} \geq \mathbf{0}, \boldsymbol{w}^{\top} \mathbf{1}=1
\end{array}$$

**Analysis**

Asboth quadratic function,$\frac{1}{2} \boldsymbol{w}^{\top} \Sigma \boldsymbol{x}$ , and afine function, $-\lambda \boldsymbol{\mu}^{\top} \boldsymbol{w} $, are convex function. And the constraint $\boldsymbol{w} \geq \mathbf{0}, \boldsymbol{w}^{\top} \mathbf{1}=1
$ are convex set. This is a convex problem.

**Code**

```matlab
sigma = [1 0.0015 -0.02;
         0.0015 1 -0.1;
         -0.02 -0.1 1];
mu = [0.001,0.05,0.005]';
n = 3;

lambda  = zeros(11,1);
for i = 1:11
    lambda(i) = (i-1)*10^4;
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
plot(expected_return,volatility)
```

**Figure**

### Problem3

**Code**

```matlab
n = randi([1,100],1);
m = randi([1,100],1);
X = randi(m,n);
y = randi(m,1);
cvx_begin
    variable x(n)
    minimize( norm( y-X*x, 2 )^2 )
    subject to
        for i = 1:n-1
            x(i) <= x(i+1)
        end
cvx_end
plot(x)
xlabel("index i",'fontsize', 16)
ylabel('Optimal solution element ${\beta_{i}}$','Interpreter','latex','fontsize', 16)
```

For $n=73, m=44$, the optimal solution $\beta^*$ is

```matlab
beta =

    0.0003
    0.0018
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0065
    0.0130
    0.0130
    0.0130
    0.0130
```

the optimal value is +42.0211