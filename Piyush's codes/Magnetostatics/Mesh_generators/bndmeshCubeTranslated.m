function mesh = bndmeshCubeTranslated(N,L,T)
    m = mshCube(N,L);
    mesh = m.bnd;
    mesh = mesh.translate(T);
end