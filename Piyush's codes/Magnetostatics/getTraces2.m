% function that returns the two traces of a vector potential 
% A = [z; x; y];
% divA = 0
% \Delta A = 0
% curlA = [1;1;1];

function [TDA,TNA,curlAdotn] = getTraces2(Gamma)
    [X,~] = Gamma.qud;
    
    % Computing the fields on the points X
    A = [X(:,3) X(:,1) X(:,2)];
    curlA = ones(size(X));
    
    % Taking traces of the fields at the points X
    normals = Gamma.qudNrm;
    TDA = A - dot(A,normals,2).*normals;
    TNA = cross(curlA,normals,2);
    curlAdotn = dot(curlA,normals,2);
end