function [p] = invPerm(I)
% Given a list I containing a permutation sigma of the numbers from 1 to n, find
% the list p containing the permutation phi 
% sigma o phi(i) = i


p = 0*I;
for i = 1:length(I)
    p(I(i)) = i;
end


end

