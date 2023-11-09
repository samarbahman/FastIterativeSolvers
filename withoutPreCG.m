%% FIS Priject 1 - Samar Bahman 416255
%% Clear Space
clc;
clear;
close all;
%% CG Algorithm

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
r0Norm = norm(r);
p = r;

% Start recording the run time
tic;

% Algorithm
tol = 1e-8;
criteria = tol + 1;
itnum = 1;
while criteria > tol
    rho = dot(r, r);
    Ap = mcsr(data,p);
    alph = rho/dot(Ap, p);
    x = x + alph*p;
    e = x - x_known;
    Ae = mcsr(data,e);
    normeA(itnum) = sqrt(dot(Ae, e));
    r = r - alph*Ap;
    normr(itnum) = norm(r);
    criteria = norm(r)/r0Norm;
    bet = dot(r, r)/rho;
    p = r + bet*p;
    itnum = itnum + 1;
end

% Finish recording the runtime
toc;

% Plot Results
i = 1:1:itnum-1;
slg = semilogy(i,normeA,i,normr);
slg(1).LineWidth = 1;
slg(1).Color = "magenta";
slg(2).Color = "cyan";
slg(2).LineStyle = ":";
grid on
xlabel('Number of Iterations','Interpreter','latex')
legend('$||e||_A$','$||r||_2$','Interpreter','latex')
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