% Function to find row wise intersection and diff
function [intersection,diffA,diffB] = rowWiseIntersectionDiff(A,B)
    % Convert to cell first to apply cellfun
    cellA = mat2cell(A,ones(size(A,1),1));
    cellB = mat2cell(B,ones(size(B,1),1));

    % Apply intersection element wise
    intersection = cellfun(@intersect,cellA,cellB,'UniformOutput',false);

    % Apply diff element wise
    diffA = cellfun(@setdiff,cellA,intersection,'UniformOutput',false);
    diffB = cellfun(@setdiff,cellB,intersection,'UniformOutput',false);

    % Converting back to matrix
    intersection = cell2mat(intersection);
    diffA = cell2mat(diffA);
    diffB = cell2mat(diffB);

end