function [nxgradU,Wh] = nxgrad_P1(Vh)

[gradU,Xh] = grad_P1(Vh);
[ncross,Wh] = nx_RWGNED(Xh);

nxgradU = ncross*gradU;


end