function [B] = repeatColumns(A,N)

% Creates a matrix B such that for all n = 1..N, B(:,n:N:end) = A
% Thus, B contains each column of A repeated N times

s = size(A);
B = reshape(repmat(A,N,1),s(1),N*s(2));


end

