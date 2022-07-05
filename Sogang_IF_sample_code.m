% Problem Set 1: Sample Code
% Yun Jung Kim

% Parameter values
beta  = 0.95;   % discount factor
delta = 0.1;   % depreciation rate
alpha = 0.3;    % capital share
rho   = 0.95;   % AR1 coefficient of technology
omega = 2; %relative risk aversion
friela = 2; %Friscsh elasticity of labor supply

% Variable types
n = 1;          % Numer of controls/free variables
m = 2;          % Number of states/exogenous
x = 0;          % Number of redundant variables (i.e. labor supply)

% Steady state values
astar = 1;
rstar = 1/beta - (1 - delta);
kstar = (1/(alpha*beta)-(1-delta))^(1/(alpha-1));
cstar = kstar^alpha - delta*kstar;

cpos =1;
kpos =2;
apos =3;

B1 = zeros(n + m, n + m);
B2 = zeros(n + m, n + m);

B1(1,cpos) = -(omega*cstar)/(cstar-(1/friela));
B1(1,kpos) = beta*rstar*(alpha -1);
B1(1,apos) = beta*rstar;
B1(2,kpos) = 1;
B1(3,apos) = 1;

B2(1,cpos) = -(omega*cstar)/(cstar-(1/friela));
B2(2,cpos) = -cstar/kstar;
B2(2,kpos) = 1/beta;
B2(2,apos) = kstar^(alpha-1);
B2(3,apos) = rho;

M = inv(B1)*B2;
H11 = M(1:n,1:n);
H12 = M(1:n,n+1:n+m);
H21 = M(n+1:n+m,1:n);
H22 = M(n+1:n+m,n+1:n+m);

% Compute eigevalues and eignevectors
[L Gamm j] = eig_order(M);
invGamm = inv(Gamm);
G11 = invGamm(1:j-1,1:n);
G12 = invGamm(1:j-1,n+1:n+m);
G21 = invGamm(j:n+m,1:n);
G22 = invGamm(j:n+m,n+1:n+m);

% Policy function
P = -inv(G21)*G22;

%VAR form
A = [zeros(n,n), -inv(G21)*(G22)*(-H21*inv(G21)*G22+H22); zeros(m,n), -H21*inv(G21)*G22+H22];
B = [zeros(n,n), -inv(G21)*(G22); zeros(m,n),eye(m,m)];
delta = [0;0;1];

% Calculate impulse responses to a 1 percent shock to technology
T = 200; % Number of periods
IRF=zeros(T,n+m);
IRF(1,:) = (B*delta)';
for t=2:T
    IRF(t,:) = (A*IRF(t-1,:)')';
end

%% Figures
set(0, 'DefaultAxesFontSize', 15);
set(0, 'DefaultTextFontSize', 15);
set(0, 'DefaultLineLinewidth', 2);

figure(1)
plot(IRF(:,1))
title('Consumption')
ylabel('Percent Deviation from Steady State')
xlabel('Periods from Technology Shock')

figure(2)
plot(IRF(:,2))
title('Capital')
ylabel('Percent Deviation from Steady State')
xlabel('Periods from Technology Shock')

figure(3)
plot(IRF(:,3))
title('Technology')
ylabel('Percent Deviation from Steady State')
xlabel('Periods from Technology Shock')