%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Step4.m
%
% Curved cantilever beam subject to a tip load using Timoshenko frame elements.
%
% (c) 2015 MAE M168
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; close all; clc;

%% Input Parameters
NEL = 20; % number of elements
Radius = 10;

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
L = 1;

% Calculations of material arrays
kGA = k.*G.*A;
EI = E.*I;
EA = E.*A;

% Loads
P = -10^-4;
qAxial = zeros(NEL,1);
qTransverse = zeros(NEL,1);

%% Generate required arrays for assembly function. 
% Within this section the assembly function is called once to generate a single
% plot of the initial and deformed configurations.
NNodes = NEL+1;
EQN = zeros(3,NNodes);
CNX = zeros(2,NEL);

nActiveDoF = 3*(NNodes - 1);

EQN(4:end) = 1:nActiveDoF;
CNX(1,:) = 1:NNodes-1;
CNX(2,:) = 2:NNodes;

% Generate nodal positions
theta = linspace(pi,pi/2,NNodes);
x = Radius*cos(theta)+Radius;
z = Radius*sin(theta);
X = [x; z];
D1 = zeros(size(EQN));

% Call assembly function to compute W, R, K
[W, R, K] = TimoshenkoAssembly(EA,EI,kGA,CNX,EQN,X,D1,qAxial,qTransverse);

Q = zeros(nActiveDoF,1);
Q(end-1) = P;

D = K\(Q-R);

% New nodal positions are computed by adding displacements D
xNew = x' + ([0;D(1:3:end)]);
zNew = z' + ([0;D(2:3:end)]);

% Plot initial and deformed configurations
figure(1)
plot(x,z,'-s',xNew,zNew,'-o','linewidth',1.25)
set(gca,'fontsize',12)
legendHandle = legend('Initial','Deformed');
set(legendHandle,'interpreter','latex','fontsize',16,'edgecolor','w',...
    'location','southeast')
xlabel('$x$','interpreter','latex','fontsize',20);
ylabel('$z$','interpreter','latex','fontsize',20);
box off
axis equal

%% Generate data for tip deflection and rotation vs log of element size
nPts = 20; 
maxElementExponent = 3;
NELSpace = round(logspace(0,maxElementExponent,nPts));

wTimo = zeros(nPts,1);
thetaTimo = zeros(nPts,1);
count = 0;

for i = NELSpace
    
    % All parameters are recomputed as NEL changes each iteration
    NEL = i;
    NNodes = i + 1;
    EQN = zeros(3,NNodes);
    CNX = zeros(2,i);
    
    nActiveDoF = 3*(NNodes - 1);
    
    EQN(4:end) = 1:nActiveDoF;
    CNX(1,:) = 1:NNodes-1;
    CNX(2,:) = 2:NNodes;
    
    % Cross sectional properties
    I = b*h^3/12*ones(NEL,1);
    A = b*h*ones(NEL,1);
    
    % Material Properties
    E = 1*ones(NEL,1);
    G = E/(2*(1+nu));
    
    % Calculations of material arrays
    kGA = k.*G.*A;
    EI = E.*I;
    EA = E.*A;
    
    % Loads
    qAxial = zeros(NEL,1);
    qTransverse = zeros(NEL,1);

    % Nodal positions
    theta = linspace(pi,pi/2,NNodes);
    x = Radius*cos(theta)+Radius;
    z = Radius*sin(theta);
    X = [x; z];
    
    D1 = zeros(size(EQN));
    
    [~, R, K] = TimoshenkoAssembly(EA,EI,kGA,CNX,EQN,X,D1,qAxial,qTransverse);
    
    % Solve for unknown displacements
    Q = zeros(nActiveDoF,1);
    Q(end-1) = P;
    D = K\(Q-R);
    
    % Add tip deflection and rotation to arrays
    count = count + 1;
    wTimo(count) = D(end-1);
    thetaTimo(count) = D(end);

end

elementSize = 1./NELSpace; % equivalent to elementLength/totalLength

% Generate plot showing convergence 
figure(2)
semilogx(elementSize,wTimo,'s--',elementSize,thetaTimo,'o-','linewidth',1.25)
set(gca,'fontsize',12)
legendHandle = legend('$w$','$\theta$');
set(legendHandle,'interpreter','latex','fontsize',16,'edgecolor','w',...
    'location','best')
xlabel('$\log{\left(\frac{\ell}{L}\right)}$','interpreter','latex',...
    'fontsize',20);
ylabel('Tip Deformation','interpreter','latex','fontsize',20);
box off






