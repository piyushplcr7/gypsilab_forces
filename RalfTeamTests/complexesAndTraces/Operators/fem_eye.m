function[I] = fem_eye(fe)

I = spdiags(ones(fe.ndof,1),0,fe.ndof,fe.ndof);

end