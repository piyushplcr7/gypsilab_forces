function [C] = alternateLines(A,B)

% See alternateColumns.

C = alternateColumns(A.',B.').';


end

