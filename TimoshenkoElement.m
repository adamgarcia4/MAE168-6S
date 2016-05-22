function [w, r, k] = TimoshenkoElement(EA,EI,kGA,x,d,qu,qw)

% Length and angle properties
dx = x(3) - x(1);
dy = x(4) - x(2);

L = sqrt(dx^2 + dy^2);

c = dx/L;
s = dy/L;

% Transformation Matrix
T1 = [c -s 0; s c 0; 0 0 1];
T = [T1 zeros(3); zeros(3) T1];

% Compute stiffness matrix
kLocal = [EA/L, 0, 0, -EA/L, 0, 0;
    0, kGA/L, kGA/2, 0, -kGA/L, kGA/2;
    0, kGA/2, kGA*L/3 + EI/L, 0, -kGA/2, kGA*L/6 - EI/L;
    -EA/L, 0, 0, EA/L, 0, 0;
    0, -kGA/L, -kGA/2, 0, kGA/L, -kGA/2;
    0, kGA/2, kGA*L/6 - EI/L, 0, -kGA/2, kGA*L/3 + EI/L];

k = T*kLocal*T';

% Distributed load vector
fqLocal = [qu*L/2 qw*L/2 0 qu*L/2 qw*L/2 0]';
fq = T*fqLocal;

% Residual force
r = k*d - fq;

% Energy
w = .5*d'*k*d - fq'*d;

end
