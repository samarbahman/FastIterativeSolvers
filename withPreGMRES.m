%% FIS Priject 1 - Samar Bahman 416255
%% Clear Space
clc;
clear;
close all;
%% GMRES Algorithm (with preconditioning)
% Read the MCSR matrix data from a text file
T = readtable('gmres_test_msr.txt');
T = table2cell(T);
data = cell2mat(T);
JM = data(2:end,1);
VM = data(2:end,2);
n = data(1,1);              % matrix size
x_known = ones(n,1);        % desired x
b = mcsr(data,x_known);     % desired b

% Initialization
x0 = zeros(n,1);
m = 600; 
tol = 1e-8;

% Preconditioned GMRES Method
preCond = input('Enter 0 for none, 1 for Jacobi, and 2 for Guass-Seidel: ');

M = zeros(n);

switch preCond
    case 0                              % No Preconditioner
        for i = 1:n
            M(i,i) = 1;
        end
%         [x, e] = gmRES(data, M, b, x0, m, tol);
%         save("e_none.mat","e")
    case 1                              % Jacobi Preconditioner
        for i = 1:n
            M(i,i) = VM(i);
        end
%         [x, e] = gmRES(data, M, b, x0, m, tol);
%         save("e_Jacobi.mat","e")
    case 2                              % Gauss-Seidel Preconditioner
        for i = 1:n
            M(i,i) = VM(i);
            i1 = JM(i);
            i2 = JM(i+1) - 1;
            for j = i1:i2
                col = JM(j);
                if col < i
                    M(i,col) = VM(j);
                else
                    break;
                end
            end
        end
%         [x, e] = gmRES(data, M, b, x0, m, tol);
%         save("e_GS.mat","e")
end

% Restarted GMRES Method
reStart = input('Enter 1 for full and 2 for restarted: ');

% Start recording the run time
tic;

switch reStart
    case 1                                  % Full GMRES Method
        [x, e] = gmRES(data, M, b, x0, m, tol);
    case 2                                  % Restarted GMRES Method
        maxIter = 1000;
        for iter = 1:maxIter
            [x, e] = gmRES(data, M, b, x0, m, tol);
            disp(iter)
            disp(e)
            x0 = x;
            if e < tol
                disp('Converged!');
                break;
            end
        end
end

% Finish recording the runtime
toc;
%% GMRES Algorithm Function
function [x, e] = gmRES(data, M, b, x, m, tol)
r0 = b - mcsr(data,x);
res = backSub(M, b);
r = backSub(M, r0);
v(:,1) = r/l2Norm(r);

% Initialization
sn = zeros(m, 1);
cs = zeros(m, 1);
e1 = zeros(m+1, 1);
e1(1) = 1;
g = l2Norm(r)* e1;
    
for k = 1:m

    % Arnoldi Algorithm
    [h(1:k+1, k), v(:, k+1)] = arnoldi(data, M, v, k);
    
    % QR Factorization using Given's Rotations
    [h(1:k+1, k), cs(k), sn(k)] = qrFact(h(1:k+1,k), cs, sn, k);
    
    % Updating the Vector
    g(k + 1) = -sn(k) * g(k);
    g(k) = cs(k) * g(k);
    error = abs(g(k + 1)) / l2Norm(res);
    e(k) = error;
    
    % Tolerance Condition
    if (error <= tol)
      break;
    end
end

% Final Result
y = backSub(h(1:k, 1:k), g(1:k));
x = x + v(:, 1:k) * y;
end
%% Arnoldi Algorithm
function [h, v] = arnoldi(data, M, v, k)
z = mcsr(data,v(:,k));
w = backSub(M, z);

for i = 1:k     
    h(i) = w' * v(:, i);
    w = w - h(i) * v(:, i);
end
h(k + 1) = l2Norm(w);
v = w / h(k + 1);
end
%% QR Factorization Function
function [h, cs_k, sn_k] = qrFact(h, cs, sn, k)
for i = 1:k-1
    temp = cs(i) * h(i) + sn(i) * h(i + 1);
    h(i+1) = -sn(i) * h(i) + cs(i) * h(i + 1);
    h(i) = temp;
end

[cs_k, sn_k] = givensRotation(h(k), h(k + 1));

h(k) = cs_k * h(k) + sn_k * h(k + 1);
h(k + 1) = 0.0;
end
%% Given's Rotation Function 
function [cs, sn] = givensRotation(v1, v2)
cs = v1 / sqrt(v1^2 + v2^2);  
sn = v2 / sqrt(v1^2 + v2^2);
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
    end
end
end
%% Dot Product Function
function dotProduct = dotProduct(vec1,vec2)
for i = 1:numel(vec1)
    product(i) = vec1(i)*vec2(i);
end
dotProduct = sum(product);
end
%% Euclidean Norm Function 
function l2Norm = l2Norm(vec)
l2Norm = sqrt(dotProduct(vec,vec));
end
%% Backwards Substitution Function
function x = backSub(U, y)
n = size(U, 1);
x = zeros(n, 1);

% Backward substitution
for i = n:-1:1
    x(i) = (y(i) - U(i, i+1:n) * x(i+1:n)) / U(i, i);
end
end