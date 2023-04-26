function mesh = meshSymTetra()
    v1 = [1 0 -1/sqrt(2)];
    v2 = [-1 0 -1/sqrt(2)];
    v3 = [0 1 1/sqrt(2)];
    v4 = [0 -1 1/sqrt(2)];

    vtcs = [v1;v2;v3;v4];
    elts = [1 2 3; 2 4 3; 3 4 1; 2 1 4];

    mesh = msh(vtcs,elts);

end