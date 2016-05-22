function [w, r, k] = BeamElement(EA,EI, x, d, q)

% Length and angle properties
dx = x(3) - x(1);
dy = x(4) - x(2);

l = sqrt(dx^2 + dy^2);

c = dx/l;
s = dy/l;

% Transformation Matrix
T1 = [c -s 0; s c 0; 0 0 1];
T = [T1 zeros(3); zeros(3) T1];

% Generate stiffness matrix
term1 = EA/l;
term2 = EI/l^3;
term3 = EI/l^2;
term4 = EI/l;

% Each row of k
k1 = [term1 0 0 -term1 0 0];
k2 = [0 12*term2 6*term3 0 -12*term2 6*term3];
k3 = [0 6*term3 4*term4 0 -6*term3 2*term4];
k4 = -k1;
k5 = -k2;
k6 = [0 6*term3 2*term4 0 -6*term3 4*term4];

% Concatenate element stiffness rows
kPrime = [k1;k2;k3;k4;k5;k6];
k = T*kPrime*T';

% Distributed load vector assuming no axial distributed load
fqPrime = [0 q*l/2 q*l^2/12 0 q*l/2 -q*l^2/12]';
fq = T*fqPrime;

% Residual Force
r = k*d - fq;

% Energy
w = .5*d'*k*d - fq'*d;

end
