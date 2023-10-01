%% FIS Priject 1 - Samar Bahman 416255
%% Clear Space
clc;
clear;
close all;
%% CG Algorithm

% Start recording the run time
tic;

% Read the MCSR matrix data from a text file
T = readtable('cg_test_msr.txt');
T = table2cell(T);
data = cell2mat(T);

n = data(1,1);              % matrix size
x_known = ones(n,1);        % desired x
b = mcsr(data,x_known);     % desired b

% Initialization
x = zeros(n,1);
r = b - mcsr(data,x);
r0Norm = l2Norm(r);
p = r;

% Algorithm
for i = 1:100000
    rho = dotProduct(r, r);
    Ap = mcsr(data,p);
    alph = rho/dotProduct(Ap, p);
    x = x + alph*p;
    %x_final(:,i) = x;
    %e(:,i) = x - x_known;
    e = x - x_known;
    Ae(:,i) = mcsr(data,e);
    normeA(i) = sqrt(dotProduct(Ae, e));
    r = r - alph*Ap;
    normr(i) = sqrt(dotProduct(r, r));
    if l2Norm(r)/r0Norm < 1e-8
        break;
    end
    bet = dotProduct(r, r)/rho;
    p = r + bet*p;
end

% Save data

save("norm.mat","normeA","normr")

% Finish recording the runtime
toc;
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