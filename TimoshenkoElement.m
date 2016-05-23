function [w, r, k] = TimoshenkoElement(EA,EI,kGA,x,d,qu,qw,Q)

if 1==1
        % Length and angle properties
    dx = x(3) - x(1);
    dy = x(4) - x(2);

    L = sqrt(dx^2 + dy^2);

    c = dx/L;
    s = dy/L;
    kLocal = zeros(6,6);

    % Transformation Matrix
    T1 = [c -s 0; s c 0; 0 0 1];
    T = [T1 zeros(3); zeros(3) T1];

    %Gauss Abscissa
    xi = zeros(Q,1);
    %Gauss Weights
    w = zeros(Q,1);

    %Set Abscissa and weights for given n-point gauss quadrature
    %Allows for one-point or two-point rule.
    if Q==1
       xi(1,1) = 0;
       w(1,1) = 2;
    elseif Q==2
        xi(1,1) = -0.5773502691896257;
        xi(2,1) = 0.5773502691896257;
        w(1,1) = 1;
        w(2,1) = 1;
    end


    for p=1:Q %Quadrature loop
       n1 = 0.5*(1-xi(p,1));
       n2 = 0.5*(1+xi(p,1));
       dn1dx = -1/(L);
       dn2dx = 1/(L);

       upr = d(1)*dn1dx + d(4)*dn2dx;
       vpr = d(2)*dn1dx + d(5)*dn2dx;
       th = d(3)*n1 + d(6)*n2;
       thpr = d(3)*dn1dx + d(6)*dn2dx;

       J = 0.5*L;
       %J = 2/(x(3)-x(1)); %Jacobian
       B1 = [dn1dx 0 0;0 0 dn1dx;0 -dn1dx n1];
       B2 = [dn2dx 0 0;0 0 dn2dx;0 -dn2dx n2];

       B = [B1 B2];
       D = [EA 0 0; 0 EI 0; 0 0 kGA];

       %Stiffness Matrix Contribution
       kLocal = kLocal + B'*D*B*J* w(p,1);
       %Element energy Contribution
       w = w + (((0.5* EI*thpr^2)+(0.5*kGA*(vpr-th))+(0.5*EA*(upr)^2))*J)*w(p,1);

    end

    k = T*kLocal*T';
    %w = .5*d'*k*d;
    r = k*d;
else

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




% k = T*kLocal*T';
% 
% % Distributed load vector
% fqLocal = [qu*L/2 qw*L/2 0 qu*L/2 qw*L/2 0]';
% fq = T*fqLocal;
% 
% % Residual force
% r = k*d - fq;
% 
% % Energy
% w = .5*d'*k*d - fq'*d;
% 
% end
% 


% 
% 
% 
% 
% % Compute stiffness matrix
% kLocal = [EA/L, 0, 0, -EA/L, 0, 0;
%     0, kGA/L, kGA/2, 0, -kGA/L, kGA/2;
%     0, kGA/2, kGA*L/3 + EI/L, 0, -kGA/2, kGA*L/6 - EI/L;
%     -EA/L, 0, 0, EA/L, 0, 0;
%     0, -kGA/L, -kGA/2, 0, kGA/L, -kGA/2;
%     0, kGA/2, kGA*L/6 - EI/L, 0, -kGA/2, kGA*L/3 + EI/L];
% 
% k = T*kLocal*T';
% 
% % Distributed load vector
% fqLocal = [qu*L/2 qw*L/2 0 qu*L/2 qw*L/2 0]';
% fq = T*fqLocal;
% 
% % Residual force
% r = k*d - fq;
% 
% % Energy
% w = .5*d'*k*d - fq'*d;

end


% function [w, r, k] = TimoshenkoElement(EA,EI,kGA,x,d,qu,qw,Q)
% 
% % Length and angle properties
% dx = x(3) - x(1);
% dy = x(4) - x(2);
% 
% L = sqrt(dx^2 + dy^2);
% 
% c = dx/L;
% s = dy/L;
% 
% % Transformation Matrix
% T1 = [c -s 0; s c 0; 0 0 1];
% T = [T1 zeros(3); zeros(3) T1];
% 
% % Compute stiffness matrix
% kLocal = [EA/L, 0, 0, -EA/L, 0, 0;
%     0, kGA/L, kGA/2, 0, -kGA/L, kGA/2;
%     0, kGA/2, kGA*L/3 + EI/L, 0, -kGA/2, kGA*L/6 - EI/L;
%     -EA/L, 0, 0, EA/L, 0, 0;
%     0, -kGA/L, -kGA/2, 0, kGA/L, -kGA/2;
%     0, kGA/2, kGA*L/6 - EI/L, 0, -kGA/2, kGA*L/3 + EI/L];
% 
