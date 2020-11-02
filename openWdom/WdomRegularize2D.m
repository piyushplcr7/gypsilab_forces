function Ms = WdomRegularize2D(data)

%%% INPUT ANALYSIS


if (length(data) == 4)
    X     = data{1};
    Ydom  = data{2};
    green = data{3};
    v     = data{4};
elseif length(data) == 5
    Xdom  = data{1};
    Ydom  = data{2};
    u     = data{3};
    green = data{4};
    v     = data{5};
else
    Xdom  = data{1};
    Ydom  = data{2};
    u     = data{3};
    F     = data{4};
    green = data{5};
    v     = data{6};
end


%%% INITIALIZATION
% Mesh data from Y
vtx  = Ydom.msh.vtx;
elt  = Ydom.msh.elt;
ctr  = Ydom.msh.ctr;
nrm  = Ydom.msh.nrm;
stp  = Ydom.msh.stp;
lgt  = Ydom.msh.ndv;
tau  = Ydom.msh.tgt;
nu   = cell2mat({tau,tau});
Nelt = size(elt,1);
gss = Ydom.gss;

% Quadrature data from Y
[~,Wy,elt2qud] = Ydom.qud;

% Degrees of freedom from Y
[~,elt2dof] = v.dof;
Nbas        = size(elt2dof,2);

% Quadrature data from X
if (length(data) >= 5)
    [X,Wx] = Xdom.qud;
end
Nx = size(X,1);

% Rangesearch with max(|edge|)_Y
[Ielt,Relt] = rangeSearch(X,ctr,1.5*stp(2));                  %%% DEBUG %%%
Mx          = cell(Nelt,1);

%%% RIGHT INTEGRATION WITH REGULARIZATION
for el = 1:Nelt
    % Edge data for Y
    Sel  = vtx(elt(el,:),:);
    A = Sel(1,:);
    B = Sel(2,:);
    Nel  = nrm(el,:);
    Tel  = tau(el,:);
    NUel = reshape(nu(el,:),3,2)';
    
    % Local size
    rMin = 10*lgt(el);                                       %%% DEBUG %%%
    
    % Quadratures points in interaction
    Iy = elt2qud(el,:);
    Ix = sort(Ielt{el}(Relt{el}<rMin))';
    
    % If interactions
    if ~isempty(Ix)
        %%% CORRECTION WITH SEMI-ANALYTIC INTEGRATION

        [WlogR,WrlogR,WgradlogR,sX,d,omegaX] = WdomSemiAnalyticInt2D(Ydom,X(Ix,:),A,B,el);
        
        Wq = Wy(Iy);
        sY = Ydom.r_1x(Iy);
        omegaX_1 = (1./omegaX);
        sX_sY = (sX*ones(1,gss) - ones(length(Ix),1)*sY');
        r = -sX_sY.*(omegaX_1*ones(1,gss));
        R2 = d.^2 + r.^2;
        
        logR = 1/2*log(R2);
        logR(sqrt(R2)<1e-13) = log(1e-13);
        logR(isnan(logR)) = 0;
        WlogR = WlogR - logR*Wq;     
        WlogR(isnan(WlogR)) = 0;
        
        rlogRt      = r.*logR;
        rlogRn      = d.*logR;
        rlogR1      = rlogRt*Tel(1) + rlogRn*Nel(1);
        rlogR2      = rlogRt*Tel(2) + rlogRn*Nel(2);
        rlogR3      = rlogRt*Tel(3) + rlogRn*Nel(3);
        WrlogR(:,1) = WrlogR(:,1) - rlogR1 * Wq;
        WrlogR(:,2) = WrlogR(:,2) - rlogR2 * Wq;
        WrlogR(:,3) = WrlogR(:,3) - rlogR3 * Wq;    
        WrlogR(isnan(WrlogR)) = 0;
        
        R_2           = 1./R2;
        R_2(isinf(R_2)) = 0;
        R_2(isnan(R_2)) = 0;
        gradlogRt       = r.*R_2;
        gradlogRn       = d.*R_2;
        gradlogR1       = gradlogRt*Tel(1) + gradlogRn*Tel(1);
        gradlogR2       = gradlogRt*Tel(2) + gradlogRn*Tel(2);
        gradlogR3       = gradlogRt*Tel(3) + gradlogRn*Tel(3);
        WgradlogR(:,1)  = WgradlogR(:,1) - gradlogR1 * Wq;
        WgradlogR(:,2)  = WgradlogR(:,2) - gradlogR2 * Wq;
        WgradlogR(:,3)  = WgradlogR(:,3) - gradlogR3 * Wq;
        
        V = [];
       
        
        %%% FINITE ELEMENT P0
        if strcmp(v.typ,'P0')
            % Correction
            if strcmp(green,'[log(r)]') && strcmp(v.opr,'[psi]')
                V = WlogR;
                
            elseif strcmp(green,'grady[log(r)]') && strcmp(v.opr,'n*[psi]')
                V = WgradlogR * Nel';
                
            elseif strcmp(green,'[log(r)]') && strcmp(v.opr,'n*[psi]')
                V{1} = WlogR .* Nel(1);
                V{2} = WlogR .* Nel(2);
                V{3} = WlogR .* Nel(3);
                
            elseif strcmp(green,'[log(r)]') && strcmp(v.opr,'nxgrad[psi]')
                V{1} = zeros(size(WlogR));
                V{2} = zeros(size(WlogR));
                V{3} = zeros(size(WlogR));
                
            else
                error('domRegularize2D.m : unavailable case')
            end
            
            %%% FINITE ELEMENT P1
        elseif strcmp(v.typ,'P1')
            % For each basis function
            for j = 1:Nbas
                % Next dof
                jp1 = mod(j,2) + 1;
                
                % Height from j
                hj = ((Sel(j,:)-Sel(jp1,:)) * NUel(j,:)');
                
                % Scalar product (x-yk).nuj/hj
                tmp = ( (X(Ix,1)-Sel(jp1,1))*NUel(j,1) + ...
                    (X(Ix,2)-Sel(jp1,2))*NUel(j,2) ) ./ hj;
                
                % Correction
                if strcmp(green,'[log(r)]') && strcmp(v.opr,'[psi]')
                    V(:,j) = WlogR.*tmp;
                    
                elseif strcmp(green,'[log(r)]') && strcmp(v.opr,'x*[psi]')
                    Vx        = WlogR.*tmp;
                    V{1}(:,j) = Vx .* X(Ix,1);
                    V{2}(:,j) = Vx .* X(Ix,2);
                    V{3}(:,j) = Vx .* X(Ix,3);
                
                elseif strcmp(green,'omega2[log(r)]') && strcmp(v.opr,'n*[psi]')
                    Vx        = WlogR.*tmp + WrlogR*NUel(j,:)'/hj;
                    V{1}(:,j) = omegaX_1.^2.*Vx .* Nel(1);
                    V{2}(:,j) = omegaX_1.^2.*Vx .* Nel(2);
                    V{3}(:,j) = omegaX_1.^2.*Vx .* Nel(3);
                    
                elseif strcmp(green,'[log(r)]') && strcmp(v.opr,'omegaDomega*[psi]')
                    % omega^2 dpsi + omega omega' psi
                    
                    dOmega      = Ydom.dw;
                    omegaDomegaX = zeros(length(Ix),1);
                    omega2tauNablaX = omegaX_1.^2./hj;
                    for k = 1:3
                        omegaDomegaX = omegaDomegaX + dOmega{k}(X(Ix,:)).*omegaX_1*Tel(k);
                    end
                    V(:,j)      = WlogR.*(omegaDomegaX.*tmp + omega2tauNablaX);
                    
                elseif strcmp(green,'omega2[log(r)]') && strcmp(v.opr,'grad[psi]')
                    V{1}(:,j) = omegaX_1.^2.*NUel(1)/hj .* WlogR;
                    V{2}(:,j) = omegaX_1.^2.*NUel(2)/hj .* WlogR;
                    V{3}(:,j) = omegaX_1.^2.*NUel(3)/hj .* WlogR;
         
                    
                    
                elseif strcmp(green,'[log(r)]') && strcmp(v.opr,'n*[psi]')
                    Vx        = WlogR.*tmp + WrlogR*NUel(j,:)'/hj;
                    V{1}(:,j) = Vx .* Nel(1);
                    V{2}(:,j) = Vx .* Nel(2);
                    V{3}(:,j) = Vx .* Nel(3);
                        
                elseif strcmp(green,'[log(r)]') && strcmp(v.opr,'nxgrad[psi]')
                    NxNUj     = cross(Nel,NUel(j,:));
                    V{1}(:,j) = NxNUj(1)/hj .* WlogR;
                    V{2}(:,j) = NxNUj(2)/hj .* WlogR;
                    V{3}(:,j) = NxNUj(3)/hj .* WlogR;
                    
                elseif strcmp(green,'grady[log(r)]') && strcmp(v.opr,'n*[psi]')
                    V(:,j) = tmp .* (WgradlogR * Nel');
                    
                else
                    error('WdomRegularize2D.m : unavailable case')
                end
            end
            
        else
            error('WdomRegularize2D.m : unavailable case')
        end
        
        % Matrix-Vector product
        I = repmat(Ix,1,Nbas);
        J = repmat(elt2dof(el,:),length(Ix),1);
        if iscell(V)
            Mx{el} = [I(:) J(:) V{1}(:) V{2}(:) V{3}(:)];
        else
            Mx{el} = [I(:) J(:) V(:)];
        end
        
    end
end


%%% LEFT INTEGRATION AND RIGHT REDUCTION
% Left integration matrix
Ndof = size(v.dof,1);
if (length(data) == 4)
    Mu = speye(Nx,Nx);
    Mw = speye(Nx,Nx);
    Ms = sparse(Nx,Ndof);
elseif length(data) == 5
    Mu = u.uqm(Xdom);
    Mw = spdiags(Wx,0,length(Wx),length(Wx));
    Ms = sparse(length(u),Ndof);
else
    Mu = u.uqm(Xdom);
    Mw = spdiags(F(X).*Wx,0,length(Wx),length(Wx));
    Ms = sparse(length(u),Ndof);
end

% Regularization matrix
Mx = double(cell2mat(Mx));
if isempty(Mx)
    if iscell(Mu)
        Mx = [1 1 zeros(1,length(Mu))];
    else
        Mx = [1 1 0];
    end
end

% Left integration
if iscell(Mu)
    for i = 1:length(Mu)
        Ms = Ms + Mu{i}' * Mw * sparse(Mx(:,1),Mx(:,2),Mx(:,2+i),Nx,Ndof);
    end
else
    Ms = Mu' * Mw * sparse(Mx(:,1),Mx(:,2),Mx(:,3),Nx,Ndof);
end

% Right reduction
[~,Mv] = v.unk;
Ms     = Ms * Mv;
end
