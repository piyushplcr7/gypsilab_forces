% evaluating row wise permutation
function out = evalRowWisePerm(original_matrix,permmat)
    N = size(original_matrix,1);
    d = size(original_matrix,2);

    out = arrayfun(@(row) original_matrix(row,permmat(row,:)),1:N,'UniformOutput',false);
    % converting back to matrix form

    out = reshape(cell2mat(out),[d N])';

end