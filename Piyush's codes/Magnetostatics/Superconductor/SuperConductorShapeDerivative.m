% Superconductor shape derivative for Neumann Trace input lying in RWG
% Space

function val = SuperConductorShapeDerivative(bndmesh,TnA,Vel,DVel,omega_src,J)
    
    RWG = fem(bndmesh,'RWG');

    % Kernel, z:= y-x
    kernel = @(x,y,z) sum(z.*(Vel(x) - Vel(y)), 2)./(vecnorm(z,2,2).^3)/ (4*pi);

    Nelt = bndmesh.nelt;

    [ii,jj] = meshgrid(1:Nelt,1:Nelt);

    t1mat = panel_assembly(bndmesh,kernel,RWG,RWG,ii(:),jj(:));
    % 1st term
    val = 0.5 * dot(TnA,t1mat*TnA);

    % 2nd term 
    % Defining the kernel for single layer BIO
    % KV = @(x,y,z) 1/norm(z)/4./pi;
    KV = @(x,y,z) sqrt(1./ sum(z.^2 ,2) ) /4./pi;

    DVelRWG = RWG;
    DVelRWG.opr = 'Dvel[psi]';

    t2mat = panel_assembly_shape_derivative(bndmesh,KV,DVelRWG,RWG,ii(:),jj(:),Vel,DVel);
    val = val + dot(TnA,t2mat*TnA);

    % 3rd term
    t3 = SuperConductorShapeDerivativeT3(bndmesh,TnA,Vel,omega_src,J);

    % 4th term
    t4 = SuperConductorShapeDerivativeT4(bndmesh,TnA,DVel,omega_src,J);

    val = val - t3-t4;
end