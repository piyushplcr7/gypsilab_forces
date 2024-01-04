function sd = rewrittenSDBEM(Tdu,Tnu,Gamma,omega_src,mu,mu0,Vel,DVel,J)
    jmu = mu0-mu;
    bndmesh = Gamma.msh;
    P0 = fem(bndmesh,'P0');
    P1 = fem(bndmesh,'P1');

    V = single_layer(Gamma,P0,P0);
    K = double_layer_laplace(Gamma,P0,P1);

    % Computing DHJ columns
    [X_bndmesh,W_bndmesh] = Gamma.qud;
    [Y_src,Wsrc] = omega_src.qud;
    NX = size(X_bndmesh,1);
    NY = size(Y_src,1);
    YY = repmat(Y_src,NX,1);
    XX = repelem(X_bndmesh,NY,1);
    JYY = J(YY(:,1),YY(:,2),YY(:,3));
    gradxG = 1/4/pi * (YY-XX)./vecnorm(XX-YY,2,2).^3;
    integrand = cross(gradxG,JYY,2);
    HJ1 = sum(Wsrc.*reshape(integrand(:,1),NY,NX),1)';
    HJ2 = sum(Wsrc.*reshape(integrand(:,2),NY,NX),1)';
    HJ3 = sum(Wsrc.*reshape(integrand(:,3),NY,NX),1)';

    XXMYY = XX-YY;

    DHJ1_integrand = 3/4/pi * XXMYY.*XXMYY(:,1)./vecnorm(XX-YY,2,2).^5 -1/4/pi./vecnorm(XXMYY,2,2).^3 * [1 0 0];
    DHJ2_integrand = 3/4/pi * XXMYY.*XXMYY(:,2)./vecnorm(XX-YY,2,2).^5 -1/4/pi./vecnorm(XXMYY,2,2).^3 * [0 1 0];
    DHJ3_integrand = 3/4/pi * XXMYY.*XXMYY(:,3)./vecnorm(XX-YY,2,2).^5 -1/4/pi./vecnorm(XXMYY,2,2).^3 * [0 0 1];

    DHJ1_integrand = cross(DHJ1_integrand,JYY,2);
    DHJ1_integrand = reshape(DHJ1_integrand,NY,3*NX);
    DHJ1_integral = sum(Wsrc.*DHJ1_integrand,1);
    DHJ1_integral = reshape(DHJ1_integral,NX,3);

    DHJ2_integrand = cross(DHJ2_integrand,JYY,2);
    DHJ2_integrand = reshape(DHJ2_integrand,NY,3*NX);
    DHJ2_integral = sum(Wsrc.*DHJ2_integrand,1);
    DHJ2_integral = reshape(DHJ2_integral,NX,3);

    DHJ3_integrand = cross(DHJ3_integrand,JYY,2);
    DHJ3_integrand = reshape(DHJ3_integrand,NY,3*NX);
    DHJ3_integral = sum(Wsrc.*DHJ3_integrand,1);
    DHJ3_integral = reshape(DHJ3_integral,NX,3);

    DHJrow1 = [DHJ1_integral(:,1) DHJ2_integral(:,1) DHJ3_integral(:,1)];
    DHJrow2 = [DHJ1_integral(:,2) DHJ2_integral(:,2) DHJ3_integral(:,2)];
    DHJrow3 = [DHJ1_integral(:,3) DHJ2_integral(:,3) DHJ3_integral(:,3)];

    Vel_X_bndmesh = Vel(X_bndmesh);
    DHJvel = [dot(DHJrow1,Vel_X_bndmesh,2) dot(DHJrow2,Vel_X_bndmesh,2) dot(DHJrow3,Vel_X_bndmesh,2)];
    normals = Gamma.qudNrm;
    HJ = [HJ1 HJ2 HJ3];
    HJdotnatGamma = dot(HJ,normals,2);
    HJdotncoeffs = proj(HJdotnatGamma,Gamma,P0);

    DHJveldotn = dot(DHJvel,normals,2);
    DHJveldotn_coeffs = proj(DHJveldotn,Gamma,P0);
    
    phi = -Tnu; v = Tdu;
    dl1ds = phi' * V * DHJveldotn_coeffs; 
    dl2ds = v' * mass_matrix(Gamma,P1,P0) * DHJveldotn_coeffs;
    
    DVel1 = DVel{1}(X_bndmesh);
    DVel2 = DVel{2}(X_bndmesh);
    DVel3 = DVel{3}(X_bndmesh);

    divVel = DVel1(:,1)+DVel2(:,2)+DVel3(:,3);

    divVel_v_coeffs = proj(divVel.*reconstruct(v,Gamma,P1),Gamma,P1);

    dl3ds = HJdotncoeffs' * K * divVel_v_coeffs + DHJveldotn_coeffs' * K * v;

    r1 = -0.5 * jmu * sum(W_bndmesh.*vecnorm(HJ,2,2).^2 .*dot(Vel_X_bndmesh,normals,2),1);
    r34 = -jmu^2/mu * (HJdotncoeffs' * V * DHJveldotn_coeffs );

    sd = -mu0*(-jmu/mu * dl1ds+ jmu/2/mu0 * dl2ds - jmu/mu0 * dl3ds) + r1 + r34;
end