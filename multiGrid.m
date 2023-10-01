%% FIS Priject 2 - Samar Bahman 416255
%% Clear Space
clc;
clear;
close all;
%% Multi Grid Algorithm

% Initializing the grid
n = 9;

N = 2^n;
h = 1/N;
x = 0:h:1;
y = 0:h:1;

N_c = 2^(n-1);
h_c = 1/N_c;
x_c = 0:h_c:1;
y_c = 0:h_c:1;

% Initializing the RHS
f = zeros(numel(x), numel(y));

for j = 1:numel(y)
    for i = 1:numel(x)
        f(i,j) = 8*pi^2*sin(2*pi*x(i))*sin(2*pi*y(j));
    end 
end

% Initialization of the parameters
u_0 = zeros(numel(x), numel(y)); 
gamma = 1;
nu1 = 4;
nu2 = 1;
r_0 = abs(max(f,[],'all'));
m = 100;

% Start recording the run time
tic;

% Algorithm
for iter = 1:m
    [u_out, r_inf] = multiGrid(u_0, f, gamma, nu1, nu2);
    u_0 = u_out;
    % Convergence criteria
    convergence(iter) = r_inf/r_0;
    if convergence(iter) < 1e-10
        break;
    end
end

% Finish recording the runtime
toc;

% Save data
% save('relativeresidual9441.mat','convergence')
%% GSLEX Function
function u = gsLEX(u, f, nu)
N = length(u);
h = 1/(N-1);
x = 0:h:1;
y = 0:h:1;
for iter = 1:nu
    for j = 2:numel(y)-1
        for i = 2:numel(x)-1
            u(i,j) = (h^2*f(i,j) + u(i-1,j) + u(i+1,j) + u(i,j-1) + u(i,j+1)) / 4;
        end
    end
end
end
%% Residual Function
function [r, r_inf] = residual(u, f)
N = length(u);
h = 1/(N-1);
x = 0:h:1;
y = 0:h:1;
r = zeros(numel(x), numel(y)); 
for j = 2:numel(y)-1
    for i = 2:numel(x)-1
        r(i,j) = f(i,j) + (u(i-1,j) - 2*u(i,j) + u(i+1,j) + u(i,j-1) - 2*u(i,j) + u(i,j+1)) / h^2;
    end
end
r_inf = abs(max(r,[],'all'));
end
%% Restriction Function
function u_2h = restrict(u)
N = (length(u) - 1)/2;
h = 1/N;
x = 0:h:1;
y = 0:h:1;
u_2h = zeros(numel(x), numel(y));
for i = 2:numel(x)-1
    ii = 2*i - 1;
    for j = 2:numel(y)-1
        jj = 2*j - 1;
        u_2h(i,j) = (u(ii-1,jj-1) + 2*u(ii,jj-1) + u(ii+1,jj-1) + 2*u(ii-1,jj) + 4*u(ii,jj) + 2*u(ii+1,jj) + u(ii-1,jj+1) + 2*u(ii,jj+1) + u(ii+1,jj+1)) / 16;
    end
end
end
%% Prolongation Function
function u_h = prolong(u_2h)
N = (length(u_2h) - 1)*2;
h = 1/N;
x = 0:h:1;
y = 0:h:1;
u_h = zeros(numel(x), numel(y));
N_c = length(u_2h) - 1;
h_c = 1/N_c;
x_c = 0:h_c:1;
y_c = 0:h_c:1;
for i = 2:numel(x_c)-1
    ii = 2*i - 1;
    for j = 2:numel(y_c)-1
        jj = 2*j - 1;
        u_h(ii-1,jj-1) = u_h(ii-1,jj-1) + u_2h(i,j)/4;
        u_h(ii,jj-1) = u_h(ii,jj-1) + u_2h(i,j)/2;
        u_h(ii+1,jj-1) = u_h(ii+1,jj-1) + u_2h(i,j)/4;
        u_h(ii-1,jj) = u_h(ii-1,jj) + u_2h(i,j)/2;
        u_h(ii,jj) = u_h(ii,jj) + u_2h(i,j);
        u_h(ii+1,jj) = u_h(ii+1,jj) + u_2h(i,j)/2;
        u_h(ii-1,jj+1) = u_h(ii-1,jj+1) + u_2h(i,j)/4;
        u_h(ii,jj+1) = u_h(ii,jj+1) + u_2h(i,j)/2;
        u_h(ii+1,jj+1) = u_h(ii+1,jj+1) + u_2h(i,j)/4;
    end
end
end
%% Multigrid Function
function [u_out, r_inf] = multiGrid(u, f, gamma, nu1, nu2)
u = gsLEX(u, f, nu1);
[r, r_inf] = residual(u, f);
res = restrict(r);
rRestrict = - res;
l = length(rRestrict);
e = zeros(l);
if l < 4
    e = gsLEX(e, rRestrict, 1);
else
    for g = 1:gamma
    e = multiGrid(e, rRestrict, g, nu1, nu2);
    end
end
e_l = prolong(e);
u_out = gsLEX(u - e_l, f, nu2);
end