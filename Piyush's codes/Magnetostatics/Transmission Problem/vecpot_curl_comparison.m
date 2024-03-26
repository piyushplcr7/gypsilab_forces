% Comparing vecpot curl computation routines
N_src = 1000;
R = 2;
r = .5;
[J,mesh_src] = get_torus_source(N_src,R,r);
omega_src = dom(mesh_src,3);

X = [2 3 4;
    5 6 7];

HJ_old = compute_vecpot_curl(J,omega_src,X);

Hj_new = computeVecpotCurlTorus(1,R,r,X);