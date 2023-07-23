% Computes the gramian matrix for different tuples of a,b,c at once
% 
% Jacobian for single vectors a,b and c (3X1 vectors)
% E = [b-a c-a] (3X2)
%
% G = E' * E; (2X2)
%
% invG = (E' * E)^-1; (2X2)
%
% DCV = E * invG;
%
% For vectorization, A,B and C are assumed of sizes (NX3) where N is the 
% no. of different tuples
function [DCVx,DCVy] = GramianVectorized(A,B,C)
    BmA = B-A;
    CmA = C-A;

    % G = E' * E = [ alpha beta; beta gamma]
    alpha = dot(BmA,BmA,2);
    beta = dot(BmA,CmA,2);
    gamma = dot(CmA,CmA,2);

    detG = alpha.*gamma-beta.^2;

    % G^-1 = [al be; be ga];
    al = gamma./detG;
    be = -beta./detG;
    ga = alpha./detG;

    % Computing DCVx and DCVy
    DCVx = al .* BmA + be .* CmA;
    DCVy = be .* BmA + ga .* CmA;
    
    

end