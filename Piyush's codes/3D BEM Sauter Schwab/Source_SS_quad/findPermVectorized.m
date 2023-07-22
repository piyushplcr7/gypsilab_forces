% Finds permutation matrix row wise
function permmat = findPermVectorized(permuted_matrix, original_matrix)
    N = size(permuted_matrix,1);
    d = size(permuted_matrix,2);
    [~,permmat] = arrayfun(@(row) ismember(permuted_matrix(row,:),original_matrix(row,:)),1:N, 'UniformOutput', false);
    % Convert to matrix format
    permmat = reshape(cell2mat(permmat),[d N])';

end