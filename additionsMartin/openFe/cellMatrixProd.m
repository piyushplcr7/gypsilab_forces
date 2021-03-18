function[C] = cellMatrixProd(A,B,op)

if ~exist('op','var')||isempty(op)
    op = @(a,b)(a.*b);
end

sA = size(A);
sB = size(B);

assert(sA(2) == sB(1),'Incompatible dimensions');

C = cell(sA(1),sB(2));

for i = 1:sA(1)
    for j = 1:sB(2)
        C{i,j} = op(A{i,1},B{1,j});
        for k = 2:sB(1)
            C{i,j} = C{i,j} + op(A{i,k},B{k,j});
        end
    end
end

end