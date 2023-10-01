%% FIS Project 3 - Samar Bahman 416255
%% Clear Space
clc;
clear;
close all;
%% Power Iteration Algorithm

% Read the MCSR matrix data from a text file
T = readtable('cg_test_msr.txt');
T = table2cell(T);
data = cell2mat(T);

% Start recording the run time
tic;

% Power Iteration
tol_power = 1e-8;
[lambd_power, e_power] = powerIt(data,tol_power);

% Finish recording the runtime
toc;

% Plot
itnum = 1:1:numel(e_power) - 1;
error = e_power(2:end);
figure;
slg = semilogy(itnum,error);
slg.LineWidth = 1;
slg.Color = "cyan";
grid on
xlabel('Number of Iterations','Interpreter','latex')
ylabel('$|\lambda_{(k)} - \lambda_{(k-1)}|$','Interpreter','latex')
%% Lanczos Algorithm
m = [30, 50, 75, 100];
tol = [1e-2, 1e-4, 1e-6, 1e-10];

for i = 1:numel(m)

    % Start recording the run time
    tic;
    
    % Lanczos Algorithm
    [lambd_lanczos, e_lanczos] = lanczos(data,m(i),tol(i));
    
    % Finish recording the runtime
    toc;

    % Save Data
    save(['data_' num2str(i) '.mat'],'lambd_lanczos', 'e_lanczos')

    % Plot
    itnum = 1:1:numel(e_lanczos) - 1;
    error = e_lanczos(2:end);
    figure;
    slg = semilogy(itnum,error);
    slg.LineWidth = 1;
    slg.Color = "magenta";
    grid on
    grid minor
    xlabel('Number of Iterations','Interpreter','latex')
    ylabel('$|\lambda_{(k)} - \lambda_{(k-1)}|$','Interpreter','latex')
    title(['$ m = $' num2str(m(i))],'Interpreter','latex')
    hold off

end
%% MCSR Storage Format Function
function y = mcsr(data, x)
n = data(1,1);
JM = data(2:end,1);
VM = data(2:end,2);

% Initialize the results array
y = zeros(n,1);

for i = 1:n
    y(i) = VM(i)*x(i);                      % Calculating the diagonal elements
    i1 = JM(i);
    i2 = JM(i+1) - 1;
    for j = i1:i2
        y(i) = y(i) + VM(j) * x(JM(j));     % Calculating the lower half with CRS format
        y(JM(j)) = y(JM(j)) + VM(j) * x(i); % Calculating the upper half with CSC format
    end
end
end
%% Power Iteration Function
function [lambda, e] = powerIt(data,tol)
e = 1;
k = 2;
n = data(1,1); 
q = ones(n,1) / sqrt(n);
while e > tol
    z = mcsr(data,q);
    q = z / l2Norm(z);
    lambda(k) = q' * mcsr(data,q);
    e(k) = abs(lambda(k) - lambda(k - 1));
    k = k + 1;
end
end
%% Lanczos Algorithm Function
function [lambd, e] = lanczos(data,k,tol)
n = data(1,1); 
v(:,1) = zeros(n,1);
v(:,2) = ones(n,1)/sqrt(n);
bet(1) = 0;
for i = 2:k + 1
    w = mcsr(data,v(:,i)) - bet(i - 1) * v(:,i - 1);  
    alph(i) = v(:,i)' * w;
    w = w - alph(i) * v(:,i);
    bet(i) = l2Norm(w);
    v(:,i + 1) = w / bet(i);
end
T = zeros(k,k);
T(1:1 + k:k * k) = alph(2:k + 1);
T(k + 1:1 + k:k * k) = bet(2:k);
T(2:1 + k:k * k - k) = bet(2:k);

e = 1;
j = 2;
q = ones(k,1)/sqrt(k);
while e > tol
    z = T * q;
    q = z/l2Norm(z);
    lambd(j) = q' * (T * q);
    e(j) = abs(lambd(j) - lambd(j - 1));
    j = j + 1;
end
end
%% Dot Product Function
function dotProduct = dotProduct(vec1,vec2)
if size(vec1,1) == size(vec2,1)
    for i = 1:size(vec1,1)
        product(i) = vec1(i)*vec2(i);
    end
else
    disp("Dimension Mismatch")
end
dotProduct = sum(product);
end
%% Euclidean Norm Function 
function l2Norm = l2Norm(vec)
l2Norm = sqrt(dotProduct(vec,vec));
end