function out_mesh = genMeshFromScript(meshscript)
    mshh = meshscript();
    Vertices = mshh.POS;
    Elements = mshh.TRIANGLES;
    Elements = Elements(:,1:3);
    clear mshh;
    out_mesh = msh(Vertices,Elements);
end