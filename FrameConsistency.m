% Script FrameConsistency.m
%
% Check consistency of frame elements
%
% (c) 2015 MAE M168

clear; close all; clc;

% Seed Random number generator
rng('shuffle');

% Generate random positions and displacements
X = 10*(rand(4,1) - .5);
u = (rand(6,1) - .5);

% Generate random moduli and distributed loads
EA = rand;
EI = rand;
kGA = rand;
qu = rand;
qw = rand;

% Compute actual solution
[W,R,K] = TimoshenkoElement(EA,EI,kGA,X,u,qu,qw);

% approximate R and K by central-difference numerical differentiation
N = 100;
h = logspace(-12,6,N)';
Rerr = zeros(N,1);
Kerr = zeros(N,1);
Rnorm = zeros(N,1);
Knorm = zeros(N,1);

for p=1:N
    
    Rh = zeros(size(R));
    Kh = zeros(size(K));
    
    for i=1:numel(u)

        uh = u;
        % perturb forward
        uh(i) = u(i) + h(p);
        
        [Wp,Rp,~] = TimoshenkoElement(EA, EI, kGA, X, uh, qu, qw);
        Rh(i) = Wp;
        Kh(:,i) = Rp;
        
        % perturb backward
        uh(i) = u(i) - h(p);
        [Wp,Rp,~] = TimoshenkoElement(EA, EI, kGA, X, uh, qu, qw);
        Rh(i) = Rh(i) - Wp;
        Kh(:,i) = Kh(:,i) - Rp;
        
    end
    
    Rh = 0.5*Rh/h(p);
    Kh = 0.5*Kh/h(p);
    
    Rerr(p) = norm(R-Rh);
    Kerr(p) = norm(K-Kh);
    Rnorm(p) = norm(R);
    Knorm(p) = norm(K);

end

% Plot error vs. h on log-log scale
figure(1)
set(gca,'fontsize',12)
loglog(h,Rerr,h,Kerr,'linewidth',1.25);
legendHandle = legend('$\|\mathbf{r}-\mathbf{r}_h\|$',...
    '$\|\mathbf{k}-\mathbf{k}_h\|$');
set(legendHandle,'interpreter','latex','fontsize',16,'edgecolor','w')
xlabel('$h$','interpreter','latex','fontsize',20);
ylabel('$\mathrm{Error}$','interpreter','latex','fontsize',20);
box off