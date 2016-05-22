function [W, R, K] = BeamAssembly(EA, EI, CNX, EQN, X, D, q)
% Function to compute and assemble the energy, nodal residual force vector,
% and stiffness matrix of an assembly of EB frame elements.
%
% (c) 2015 MAE M168
%
% Input parameters:
% EA: (Vector, NEL x 1) Stretching modulus for the elements
% EI: (Vector, NEL x 1) Bending modulus for the elements
% CNX: (Matrix, 2 x NEL) Nodal connectivity array for the elements%
% EQN: (Matrix, 3 x NN) Global equation numbers
% X: (Matrix, 2 x NN) Nodal positions
% D: (Matrix, 3 x NN) Nodal displacements
% q: (Vector, NEL x 1) Transverse uniform load
%
% Output parameters
% W: (Scalar) Internal energy
% R: (Vector, N x 1) Internal nodal forces
% K: (Matrix, N x N) Stiffness matrix

E = size(CNX,2); % number of elements
N = sum(sum(EQN > 0)); % number of active DoF
NODEDoF = 3; % nodal DoF
NPE = 2; % nodes per element

% initialize energy, residual, stiffness
W = 0.0;
R = zeros(N,1);
K = zeros(N,N);

% loop over elements

for e=1:E
    
    % Get both element node IDs
    elementNodes = CNX(:,e);
    
    % Compute DoF locations in D, the global displacement array
    DoFLocations = [elementNodes(1)*NODEDoF-2;
        elementNodes(1)*NODEDoF-1;
        elementNodes(1)*NODEDoF;
        elementNodes(2)*NODEDoF-2;
        elementNodes(2)*NODEDoF-1;
        elementNodes(2)*NODEDoF];
    
    % Copy global displacement values into element displacement array
    d = D(DoFLocations);
    
    % Copy global positions into element position array
    x = X(([elementNodes(1)*NPE-1,elementNodes(1)*NPE,...
        elementNodes(2)*NPE-1,elementNodes(2)*NPE]));
    
    % Call BeamElement to compute element quantities
    
    [w,r,k] = BeamElement(EA(e),EI(e), x, d, q(e));
    
    % Assemble element contributions to global energy, residual, stiffness
    W = W+w;
    
    % Get global DoF corresponding to both element nodes
    globalDoF = reshape(EQN(:,elementNodes),[NODEDoF*NPE,1]);
    
    % Determine active DoF
    activeDoF = (globalDoF>0);
    
    % Remove inactive DoF so they aren't included in the assembly
    globalDoF(~activeDoF) = [];
    
    % Assembly
    R(globalDoF) = R(globalDoF) + r(activeDoF);
    K(globalDoF,globalDoF) = K(globalDoF,globalDoF) + k(activeDoF,activeDoF);
    
end

end


