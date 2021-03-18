classdef polymesh < msh % mesh of a polygonal domain
    
    properties
        polyVtx; % list of the vertices of the polygon (complex coords).
        r; % conformal mapping from the disk to the exterior of the polygon
    end
    
    methods (Static)
        function[out] = calcTheta2(vertices,density)
            % Assuming vertices are given in trigonometric order
            % In other words, polygon located at the left of the itinerary
            r =  extermap(polygon(vertices));
            prev = prevertex(r);
            Nvert = length(vertices);
            out = [];
            for i = 1:length(vertices)
                Z1 = vertices(i);
                if i==1
                    Z0 = vertices(Nvert);
                    Z2 = vertices(2);
                elseif i==Nvert
                    Z0 = vertices(Nvert-1);
                    Z2 = vertices(1);
                else
                    Z0 = vertices(i-1);
                    Z2 = vertices(i+1);
                end
                Z0 = (Z0 + Z1)/2;
                Z2 = (Z2 + Z1)/2;
                N1 = fix(density*abs(Z0 - Z1)) + 1;
                N2 = fix(density*abs(Z2 - Z1)) + 1;
                z0 = evalinv(r,Z0);
                z2 = evalinv(r,Z2);
                z1 = prev(i);
                t0 = angle(z0);
                delta1 = mod(t0 - angle(z1),2*pi);
                tht1 =  linspace(t0,t0 - delta1,N1);
                delta2 = mod(angle(z1) - angle(z2),2*pi);
                tht2 = linspace(angle(z1),angle(z1) - delta2,N2);
                out = [out,tht1, tht2];
            end
            out = unique(sort(mod(out,2*pi)));
        end
        function[out] = calcTheta(vertices,density)
            r =  extermap(polygon(vertices)); % Riemann mapping
            prev = prevertex(r);
            [args,I] = sort(mod(real(log(prev)/1i) + 4*pi,2*pi));
            vertices = vertices(I);
            Nsides =  length(vertices);
            out = [];
            for side_num = 1:(length(args)-1)
                A = vertices(side_num);
                B = vertices(side_num+1);
                l = abs(B - A);
                N = round(l*density) + 1;
                tmp = linspace(args(side_num),args(side_num + 1),N+1);
                out = [out; tmp(1:end-1).'];
            end
            A = vertices(end);
            B = vertices(1);
            l = abs(B-A);
            N = round(l*density)+1;
            tmp = linspace(args(Nsides),args(1) + 2*pi,N+1);
            out = [out; tmp(1:end-1).'];
        end
    end
    methods
        
        
        
        %% on the fly properties:
        function[out] = Nsides(this)
            out = length(this.vertices);
        end
        function[out] = prevertices(this)
            out = prevertex(this.r);
        end
        function[out] = theta(this)
            out = mod(real(-1i*log(this.z)),2*pi);
        end
        function[out] = Z(this)
            vtx = this.vtx;
            out = vtx(:,1) + 1i*vtx(:,2);
        end
        function[out] = z(this)
            out = evalinv(this.r,this.Z);
        end
        
        %% class constructor
        function[this] = polymesh(vertices,density)
            % vertices : [z1,z2,z3,...zP] vertices of the polygon given as complex
            % numbers.
            % N : number of desired mesh points per segment.
            r = extermap(polygon(vertices));
            theta = polymesh.calcTheta(vertices,density);
            z = exp(1i*theta);
            Z = r(z);
            Z = fillgaps(Z,density);
            Ntot = length(Z);
            elt = [(1:Ntot)' [(2:Ntot) 1]'];
            this@msh([real(Z(:)) imag(Z(:)) 0*real(Z(:))],elt);
            this.polyVtx = vertices;
            this.r = r;
        end
        %% Methods
        function[] = plot_sc(this)
            zz = this.z; zz = [zz(:); zz(1)]; ZZ = this.Z; ZZ = [ZZ(:); ZZ(1)];
            subplot(1,2,1);
            plot(real(ZZ),imag(ZZ),'-*');
            title('Mesh')
            axis equal;
            subplot(1,2,2);
            plot(real(zz),imag(zz),'-*');
            axis equal
            title('Preimage');
        end
        function[m] = premesh(this)
            z = this.z;
            x = real(z); y = imag(z); z= 0*x;
            m = msh([x,y,z],this.elt);
        end
    end
end