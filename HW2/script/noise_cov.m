function [] = noise_cov(x)

format short
mean_x = mode(x)
x_ = x - mean_x;
% 
% x_(:,1) = x_(:,1)*100;
% x_(:,2) = x_(:,2)*100;
% x_(:,3) = x_(:,3)*180/pi;

format long

cov = (x_'*x_)/size(x,1)
eig(cov)
subplot(3,1,1)
plot(x(:,1))
hold on 
plot(mean_x(1)*ones(1,size(x,1)))
subplot(3,1,2)
plot(x(:,2))
hold on 
plot(mean_x(2)*ones(1,size(x,1)))
subplot(3,1,3)
plot(x(:,3)*180/pi)
hold on 
plot(180*mean_x(3)*ones(1,size(x,1))/pi)
