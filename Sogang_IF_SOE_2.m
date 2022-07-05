% Problem Set 1: Sample Code
% Yun Jung Kim

% Parameter values
beta  = 0.96;   % discount factor
delta = 0.1;   % depreciation rate
alpha = 0.32;    % capital share
rho   = 0.95;   % AR1 coefficient of technology
omega = 1.445;   %  
sigma = 2;      %  
dbar = 0.7442;
psi = 0.000742;

% Variable types
n = 1;          % Numer of controls/free variables (C)
m = 3;          % Number of states/exogenous (d,K, A)
x = 2;          % Number of redundant variables (L, Y)

% Steady state values
astar = 1;
rstar = (1/beta)-1;
kstar = ((1/alpha)*((1-alpha)^((alpha-1)/(omega-1+alpha)))*((1/beta)-(1-delta)))^((omega-1+alpha)/(alpha*omega-omega-alpha+1));
lstar = ((1-alpha)^(1/(omega+alpha-1)))*(kstar^(alpha/(omega+alpha-1)));
ystar = (kstar^alpha)*(lstar^(1-alpha));
dstar = dbar;
cstar = (kstar^alpha)*(lstar^(1-alpha)) - (delta*kstar) - (rstar*dstar);

cpos =1;
dpos =2;
kpos =3;
apos =4;

b1 = zeros(n + m, n + m);
b2 = zeros(n + m, n + m);
b3 = zeros(n + m, x);
b4 = zeros(n + m, x);
F = zeros(x, n + m);

b1(1,cpos) = -(sigma*cstar)/(cstar - ((lstar^omega)/omega));
b1(1,kpos) = beta*(rstar+delta)*(alpha-1);
b1(1,apos) = beta*(rstar+delta);
b1(2,dpos) = -dstar/kstar;
b1(2,kpos) = 1;
b1(3,apos) = 1;
b1(4,cpos) = (-sigma*cstar)/(cstar - ((lstar^omega)/omega));
b1(4,dpos) = (psi*dstar)/(1+rstar);

b2(1,cpos) = (-sigma*cstar)/(cstar-((lstar^omega)/omega));
b2(2,cpos) = -(cstar/kstar);
b2(2,dpos) = -((1+rstar+psi*dstar)*dstar)/kstar;
b2(2,kpos) = 1-delta;
b2(3,apos) = rho;
b2(4,cpos) = -(sigma*cstar)/(cstar-((lstar^omega)/omega));

b3(1,1) = (sigma*(lstar^omega))/(cstar-((lstar^omega)/omega))+beta*(rstar+delta)*(1-alpha);
b3(4,1) = (sigma*(lstar^omega))/(cstar-((lstar^omega)/omega));
b4(1,1) = (sigma*(lstar^omega))/(cstar-((lstar^omega)/omega));
b4(2,2) = ystar/kstar;
b4(4,1) = (sigma*(lstar^omega))/(cstar-((lstar^omega)/omega));

F(1,kpos) = alpha/(alpha+omega-1);
F(1,apos) = 1/(alpha+omega-1);
F(2,kpos) = alpha*omega/(alpha+omega-1);
F(2,apos) = omega/(alpha+omega-1);

B1=b1+b3*F;
B2=b2+b4*F;

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
delta = [0;0;0;1];

% Calculate impulse responses to a 1 percent shock to technology
T = 200; % Number of periods
IRF=zeros(T,n+m);
IRF(1,:) = (B*delta)';
for t=2:T
    IRF(t,:) = (A*IRF(t-1,:)')';
end

X=zeros(T,x);
X(1,:) = F*IRF(1,:)';
for t=2:T
X(t,:) = F*IRF(t,:)'; 
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
plot(IRF(:,3))
title('Capital')
ylabel('Percent Deviation from Steady State')
xlabel('Periods from Technology Shock')

figure(3)
plot(IRF(:,4))
title('Technology')
ylabel('Percent Deviation from Steady State')
xlabel('Periods from Technology Shock')

figure(4)
plot(X(:,1))
title('Labor')
ylabel('Percent Deviation from Steady State')
xlabel('Periods from Technology Shock')

figure(5)
plot(X(:,2))
title('Output')
ylabel('Percent Deviation from Steady State')
xlabel('Periods from Technology Shock')

figure(6)
plot(IRF(:,2))
title('Debt')
ylabel('Percent Deviation from Steady State')
xlabel('Periods from Technology Shock')