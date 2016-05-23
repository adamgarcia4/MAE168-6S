%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step2.m
%
% Cantilever beam subject to a tip load using Timoshenko frame elements.
%
% (c) 2015 MAE M168
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% Input Parameters
NEL = 20; % number of elements
L = 10;

% Cross sectional properties
b = 1;
h = 1;
I = b*h^3/12*ones(NEL,1);
A = b*h*ones(NEL,1);

% Material Properties
nu = 0;
E = 1*ones(NEL,1);
k = 5/6;
G = E/(2*(1+nu));

% Calculations of material arrays
kGA = k.*G.*A;
EI = E.*I;
EA = E.*A;

% Loads
P = -10^-4;
qAxial = zeros(NEL,1);
qTransverse = zeros(NEL,1);

%% Generate required arrays for assembly function
NNodes = NEL+1;
EQN = zeros(3,NNodes);
CNX = zeros(2,NEL);

nActiveDoF = 3*(NNodes - 1);

EQN(4:end) = 1:nActiveDoF;
CNX(1,:) = 1:NNodes-1;
CNX(2,:) = 2:NNodes;

X = [linspace(0,L,NNodes); zeros(1,NNodes)];
D1 = zeros(size(EQN));

% Call assembly function to compute W, R, K
[W, R, K] = TimoshenkoAssembly(EA,EI,kGA,CNX,EQN,X,D1,qAxial,qTransverse,2);

% Generate concentrated force vector Q
Q = zeros(nActiveDoF,1);
Q(end-1) = P;

% Solve for unknown displacements
D = K\(Q-R);
finalD = reshape(D1,[],1) + [zeros(3,1);D];

figure(1)
plot(X(1,:),finalD(2:3:end),'s-',X(1,:),finalD(3:3:end),'o-','linewidth',1.25)
set(gca,'fontsize',12)
legendHandle = legend('$w(x)$','$\theta(x)$');
set(legendHandle,'interpreter','latex','fontsize',16,'edgecolor','w',...
    'location','southwest')
xlabel('$x$','interpreter','latex',...
    'fontsize',20);
ylabel('Deformation','interpreter','latex','fontsize',20);
box off

%% Generate Zienkiewicz/Taylor Plot From Slide 35
nPts = 15;
LengthArray = logspace(0.2,3,nPts);
wTimo = zeros(nPts,1); % Timoshenko tip deflections
wEB = zeros(nPts,1); % EB tip deflections
count = 0;

for i = LengthArray
    
    count = count + 1;
    L = i;
    X = [linspace(0,L,NNodes); zeros(1,NNodes)];
    D1 = zeros(size(EQN));
    
    [~, R, K] = TimoshenkoAssembly(EA,EI,kGA,CNX,EQN,X,D1,qAxial,qTransverse,2);
    
    D = K\(Q-R);
    wTimo(count) = D(end-1);
    
    [~, R, K] = BeamAssembly(EA, EI, CNX, EQN, X, D1, qTransverse);

    D = K\(Q-R);
    wEB(count) = D(end-1);
    
end

wTimoExact = P*LengthArray'.^3/(3*EI(1))+P*LengthArray'/kGA(1);
wEBExact = P*LengthArray'.^3/(3*EI(1));

figure(2)
semilogx(LengthArray,wTimo./wEB,'s-','linewidth',1.25)
set(gca,'fontsize',12)
xlabel('$\log{(L/h)}$','interpreter','latex','fontsize',20);
ylabel('$w_{\rm T}/w_{\rm EB}$','interpreter','latex','fontsize',20);
box off

