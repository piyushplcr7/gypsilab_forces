function [C] = alternateColumns(A,B)

% Returns C = [A(:,1), B(:,1), A(:,2), B(:,2),...,A(:,N),B(:,N)]
% where size(A,2)=size(B,2)=N. 

sA = size(A); sB = size(B);
assert(sum(abs(sA-sB))==0,'operation only possible on matrices of equal size');

C = [A;B];
C = reshape(C,sA(1),2*sA(2));



end

