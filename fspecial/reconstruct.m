function x = reconstruct(x0,num_pro,data)

x = x0; 

maxlsiter = 150;
gradToll = 1.0000e-030; % check 1 
alpha = 0.0100;
beta = 0.6000;
t0 = 1;
Itnlim = 16;	

k = 0;  % Counter
t = 1;  % Multiplication factor 

g0 = wGradient(x,num_pro,data);
dx = -g0;

while(1)
    
    f0 = objective(x,dx,0,num_pro,data);
	t = t0;
    [f1]  =  objective(x,dx,t,num_pro,data);
	lsiter = 0;
    
	while (f1 > f0 - alpha*t*abs(g0(:)'*dx(:))) & (lsiter<maxlsiter)
		lsiter = lsiter + 1;
		t = t * beta;
		[f1]  =  objective(x,dx,t,num_pro,data);
    end
    
    % lsiter
    
    if lsiter == maxlsiter
		disp('Reached max line search,.... not so good... might have a bug in operators. exiting... ');
		return;
	end

	if lsiter > 2
		t0 = t0 * beta;
    end 
    
	if lsiter<1
		t0 = t0 / beta;
    end
    
    x = (x + t*dx);
    
    % gradient calculation
	g1 = wGradient(x,num_pro,data);
	bk = g1(:)'*g1(:)/(g0(:)'*g0(:)+eps);
	g0 = g1;
	dx = - g1 + bk* dx;
	k = k + 1;
	
	if (k > Itnlim) | (norm(dx(:)) < gradToll) 
		break;
    end
end
return;

function [res] = objective(x,dx,t,num_pro,data);

x = x + t*dx;
data_x = zeros(185,num_pro);
for i = 0:(num_pro-1) 
data_x(:,i+1) = radon(x,180/num_pro*i);
end
obj = data_x - data; 
res=obj(:)'*obj(:);
% Adding up the TV function now.
dif_x = zeros(128,128);
dif_y = zeros(128,128);

for i = 2:127
    for j = 2:127
        dif_y(i,j) = x(i,j) - x(i-1,j);
        dif_x(i,j) = x(i,j) - x(i,j-1);
    end
end
tv = 0.07;
sum_x = sum(abs(dif_x(:).^2));
sum_y = sum(abs(dif_y(:).^2));
res = res + tv*(sum_x + sum_y);

% adding up the pixels outside fov
left = x(:,1:5);
right = x(:,122:128);
sum_x = sum(abs(left(:).^2));
sum_y = sum(abs(right(:).^2));
res = res + tv*(sum_x + sum_y);
rest = 0;

for i = 1:128
    for j = 1:128
        if( i+j<40 | i+j > 216 | j-i>88 | i-j > 88 | j<5 | j>122) 
            rest = rest + abs(x(i,j));
        end
    end
end
 
%{
for i = 1:39
    rest = rest + abs(x(i,40-i));
end
for i = 1:39
    rest = rest + abs(x(128 - i, 88 + i));
end
for i = 1:39
    rest = rest + abs(x(i, 88 + i));
end
for i = 1:39
rest = rest + abs(x(128 - i, 40 - i));
end
%}
res = res + tv*rest;
x = x - t*dx;
return

function [grad] = wGradient(x,num_pro,data)
fov = 1;
tv = 0.07;
data_x = zeros(185,num_pro);
for i = 0:(num_pro-1) 
data_x(:,i+1) = radon(x,180/num_pro*i);
end
obj = data_x - data;
grad = iradon(obj,180/num_pro,128);

mu_x= filter2([1 -1 0]',x);
mu_y= filter2([1 -1 0],x);
grad = grad + tv*(filter2([0 -1 1]',mu_x) + filter2([0 -1 1],mu_y));

rest_fov = zeros(128,128);

for i = 1:128
    for j = 1:128
        if( i+j<40 | i+j > 216 | j-i>88 | i-j > 88 | j<5 | j>122) 
            rest_fov(i,j) = x(i,j);
        end
    end
end

grad = grad + fov*rest_fov;
%{
for i = 1:39
    rest_x(i,40-i) = x(i,40-i) ;% - x(i,41-i);
    rest_y(i,40-i) = x(i,40-i) ;% - x(i+1,40-i);
    rest_x(128 - i, 88 + i) = x(128 - i, 88 + i);% - x(128 - i, 89 + i);
    rest_y(128 - i, 88 + i) = x(128 - i, 88 + i) ;% - x(129 - i, 88 + i);
    rest_x(i, 88 + i) = x(i, 88 + i) ;%- x(i, 89 + i);
    rest_y(i, 88 + i) = x(i, 88 + i); % - x(i+1, 88 + i);
    rest_x(128 - i, 40 - i) = x(128 - i, 40 - i);% - x(128 - i, 41 - i);
    rest_y(128 - i, 40 - i) = x(128 - i, 40 - i);% - x(129 - i, 40 - i);
end
grad = grad + tv*(rest_x + rest_y);
%}
return


%{
dif_x = zeros(128,128);
dif_y = zeros(128,128);

for i = 2:126
    for j = 2:126
        dif_y(i,j) = 2*x(i,j) - x(i-1,j)- x(i+1,j);
        dif_x(i,j) = 2*x(i,j) - x(i,j-1)- x(i,j+1);
    end
end

grad = iradon(obj,180/num_pro,128);
grad = grad + tv*(dif_x + dif_y);
%}
