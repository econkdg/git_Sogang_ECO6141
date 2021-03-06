% Problem Set 1: Sample Code

% Parameter values
beta  = 0.95;   % discount factor
delta = 0.1;   % depreciation rate
alpha = 0.3;    % capital share
rho   = 0.95;   % AR1 coefficient of technology
omega = 1.445;   %  
sigma = 2;      %  

% Variable types
n = 1;          % Numer of controls/free variables (C)
m = 2;          % Number of states/exogenous (K, A)
x = 2;          % Number of redundant variables (L, Y)

% Steady state values
astar = 1;
rstar = 1/beta - (1 - delta);
kstar = ((1 - alpha)^(1/(omega - 1)))*(alpha/rstar)^((alpha + omega - 1)/((omega - 1)*(1 - alpha)));
lstar = ((1-alpha)^(1/(omega+alpha-1)))*(kstar^(alpha/(omega+alpha-1)));
cstar = kstar^alpha - delta*kstar;
ystar = (kstar^alpha)*(lstar^(1-alpha));

cpos =1;
kpos =2;
apos =3;

b1 = zeros(n + m, n + m);
b2 = zeros(n + m, n + m);
b3 = zeros(n + m, x);
b4 = zeros(n + m, x);
F = zeros(x, n + m);

b1(1,cpos) = (sigma*cstar)/(cstar - ((lstar^omega)/omega));
b1(1,kpos) = -beta*rstar*(alpha-1);
b1(1,apos) = -beta*rstar;
b1(2,kpos) = 1;
b1(3,apos) = 1;

b2(1,cpos) = (sigma*cstar)/(cstar-((lstar^omega)/omega));
b2(2,cpos) = -cstar/kstar;
b2(2,kpos) = 1-delta;
b2(3,apos) = rho;

b3(1,1) = -(sigma*(lstar^omega))/(cstar-(lstar^omega)/omega)-beta*rstar*(1-alpha);
b4(1,1) = -(sigma*(lstar^omega))/(cstar-(lstar^omega)/omega);
b4(2,2) = ystar/kstar;

F(1,2) = alpha/(alpha+omega-1);
F(1,3) = 1/(alpha+omega-1);
F(2,2) = alpha*omega/(alpha+omega-1);
F(2,3) = omega/(alpha+omega-1);

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
delta = [0;0;1];

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
plot(IRF(:,2))
title('Capital')
ylabel('Percent Deviation from Steady State')
xlabel('Periods from Technology Shock')

figure(3)
plot(IRF(:,3))
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

