clear;clc;
cvx_begin
    variables x y
    minimize(abs(x)+abs(y))
    subject to
        y == -0.1*x +1
        
        
cvx_end

% for i = 1:11
%     cvx_begin
%         variable w(n)
%         minimize(w'*sigma*w - lambda(i)*mu'*w)
%         subject to
%             w >= 0
%             w'*ones(n,1) == 1
%     cvx_end 
%     w_optimal(i,:) = w';
% end