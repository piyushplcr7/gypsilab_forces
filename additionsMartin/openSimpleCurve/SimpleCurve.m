classdef SimpleCurve
    % Parametric curve, can be closed but must not intersect itself
    
    properties
        x; % Anonymous function of 1 argument
        y; % Anonymous function of 1 argument
        I; % Common interval of definition of x and y
        closed@logical;
        boundedSide; % Is the bounded side to the left or right of the boundary
        dx; % Anonymous function of 1 argument, (derivative of x(t))
        dy; % Anonymous function of 1 argument, (derivative of x(t))
        isnormal = false; % true only if the curve is parametrrized by arclength
        suppliedDer = false; % true if user supplied explicit derivative.
    end
    
    methods
        
        %% Constructor
        
        function[this] = SimpleCurve(xx,yy,II,bS)
            this.x = xx;
            this.y = yy;
            this.I = II;
            a = II(1); b = II(2); Ma = [xx(a);yy(a)]; Mb = [xx(b);yy(b)];
            this.closed = norm(Ma-Mb)<1e-11;
            if this.closed
                assert(logical(exist('bS','var')),'this curve is closed. You must pass the boundedSide (value "left" or "right") in argument');
                assert(ismember(bS,{'left','right'}),['This curve is closed. the value used for argument bS is incorrect. Please use one of the choices : "left" or "right". \n' ...
                    'choose left if the bounded component of the plane lies at the left of the curve, and right otherwise.'])
                this.boundedSide = bS;
            else
                this.boundedSide = 'none';
            end
            this.dx = safeDx(this.x,II);
            this.dy = safeDx(this.y,II);
        end
        
        function[this] = supplyDer(this,DX,DY)
            this.suppliedDer = true;
            this.dx = DX;
            this.dy = DY;
        end
        function[this] = forgetDer(this)
            this.dx = safeDx(this.x,this.I);
            this.dy = safeDx(this.y,this.I);
            this.suppliedDer = false;
        end
        function[b,err] = checkDer(this)
            if this.suppliedDer
                t = linspace(this.I(1),this.I(2),10000);
                DX = safeDx(this.x,this.I); DY = safeDx(this.y,this.I);
                err = norm([DX(t) - this.dx(t), DY(t) - this.dy(t)],'inf');
                b = err < 1e-3;
            else
                b = false;
                err = NaN;
            end
        end
        
        %% Display
        
        function[] = disp(this)
            if this.closed
                s1 = sprintf('Closed curve (bounded side on the %s of the boundary)',...
                    this.boundedSide);
            else
                s1 = 'Open curve';
            end
            s2 = sprintf('x,y : [%s,%s] -> R^2',...
                num2str(this.I(1)),num2str(this.I(2)));
            
            s3 = sprintf('Length : %s',num2str(this.length));
            disp(s1);
            disp(s2); %#ok
            disp(s3); %#ok
        end
        
        function[t,xt,yt] = emphPoints(this)
            L = this.length();
            s_indices = linspace(L/1e4,L-L/1e4,20);
            t = this.t_of_s(s_indices);
            xt = this.x(t);
            yt = this.y(t);
        end
        
        
        function[] = plot(this)
            L = this.length();
            s_indices = linspace(L/1e4,L-L/1e4,1e3); % Avoid endpoint pbs
            t = this.t_of_s(s_indices);
            xt = this.x(t);
            yt = this.y(t);
            plot(xt,yt,'k','LineWidth',4);
%             [~,xt,yt] = emphPoints(this);
%             hold on;
%             plot(xt,yt,'+','Markersize',9);
            axis equal
%             grid on
        end
        
        function[] = showTgt(this)
            [t,xt,yt] = emphPoints(this);
            hold on
            [~,T] = tgt(this,t(:));
            quiver(xt(:),yt(:),T(:,1),T(:,2),0.6,'r','LineWidth',2);
        end
        function[] = showNrm(this)
            [t,xt,yt] = emphPoints(this);
            hold on
            [~,N] = nrm(this,t(:));
            quiver(xt(:),yt(:),N(:,1),N(:,2),0.6,'r','LineWidth',2);
        end
        function[] = showUnitTgt(this)
            [t,xt,yt] = emphPoints(this);
            hold on
            [UT,~] = tgt(this,t(:));
            quiver(xt(:),yt(:),UT(:,1),UT(:,2),0.6,'r','LineWidth',2);
        end
        function[] = showUnitNrm(this)
            [t,xt,yt] = emphPoints(this);
            hold on
            [UN,~] = nrm(this,t(:));
            quiver(xt(:),yt(:),UN(:,1),UN(:,2),0.6,'r','LineWidth',2);
        end
        function[] = showFrenet(this)
            this.showTgt;
            this.showNrm;
        end
        function[] = showOsculatoryCircles(this)
            [t,xt,yt] = emphPoints(this);
            r = 1./this.curvature(t);
            [UT,~] = nrm(this,t);
            Nx = UT(:,1); Ny = UT(:,2);
            c = [xt(:) + Nx(:).*r(:),yt(:) + Ny(:).*r(:)];
            u = linspace(0,2*pi,100);
            for i = 1:length(t)
                hold on
                plot(c(i,1) + r(i)*cos(u),c(i,2) + r(i)*sin(u),'k--');
            end
        end
        
        function[] = plotCurvature(this)
            subplot(1,2,1)
            plot(this);
            showUnitTgt(this);
            showOsculatoryCircles(this);
            title('Osculatory circles')
            subplot(1,2,2);
            t = linspace(this.I(1),this.I(2),1000);
            plot(t,this.curvature(t));
            xlabel('t')
            ylabel('Curvature')
            title('Curvature')
        end
        
        
        %% Elementary manipulation
        function[l] = bnd(this)
            if this.closed
                l = [];
            else
                l = zeros(2,3);
                for i = 1:2
                    l(i,:) = [this.x(this.I(i)),this.y(this.I(i)),0];
                end
                
            end
            
        end
        
        
        function[c] = portion(this,J)
            Iold = this.I;
            aa = J(1); bb = J(2);
            c = this;
            c.I(1) = aa;
            c.I(2) = bb;
            if this.closed
                 gap = max(abs(Iold(1) - c.I(1)),abs(Iold(2) - c.I(2)));
                 if gap > 1e-6
                     c.closed = false;
                 else
                     c.I = Iold;
                 end
            end
        end
        
        function[newC] = affineReparam(this,II)
            c = II(1);
            d = II(2);
            assert(c<d,'first use reverse to modify travel direction');
            a = this.I(1); b = this.I(2);
            X = this.x; Y = this.y;
            
            xx = @(t)(X(a + (t-c)/(d-c)*(b-a)));
            yy = @(t)(Y(a + (t-c)/(d-c)*(b-a)));
            
            newC = SimpleCurve(xx,yy,[c,d],this.boundedSide);
            
            if this.suppliedDer
                DX = this.dx; DY = this.dy;
                newC = supplyDer(newC,...
                    @(t)(DX(a + (t-c)/(d-c)*(b-a)))*(b-a)/(d-c),...
                    @(t)(DY(a + (t-c)/(d-c)*(b-a)))*(b-a)/(d-c));
            end
        end
        function[c] = reverse(this)
            % Reverse the direction of travel
            X = this.x; Y = this.y;
            a = this.I(1);
            b = this.I(2);
            xx = @(t)(X(a+b - t));
            yy = @(t)(Y(a+b - t));
            if this.closed
                switch this.boundedSide
                    case 'left'
                        bS = 'right';
                    case 'right'
                        bS = 'left';
                end
            else
                bS = 'none';
            end
            c = SimpleCurve(xx,yy,this.I,bS);
            if this.isnormal
                c.isnormal = true;
            end
            if this.suppliedDer
                DX = this.dx; DY = this.dy;
                c = supplyDer(c,...
                    @(t)(DX(a+b - t)*(-1)),...
                    @(t)(DY(a + b - t)*(-1)));
            end
        end
        
        
        %% Tangent and normals
        
        
        function[UT,T] = tgt(this,t)
            Tx = this.dx(t);
            Ty = this.dy(t);
            T = [Tx(:), Ty(:)];
            UT = T./(sqrt(Tx(:).^2 + Ty(:).^2)*[1 1]);
        end
        function[UN,N] = nrm(this,t)
            % Convention : [T(x),N(x)] is direct, where T and N are the
            % tangent and normal vectors at x.
            [UT,T] = this.tgt(t);
            UN = [-UT(:,2),UT(:,1)];
            N = [-T(:,2),T(:,1)];
        end
        
        %% Arclength
        
        function[l] = length(this)
            if this.isnormal
                l = this.I(2);
            else
                l = integral(@(t)(sqrt(this.dx(t).^2 + this.dy(t).^2)),...
                    this.I(1),this.I(2));
            end
        end
        function[Gs] = s(this,t)
            filled = false;
            % Compute arclength at query points t.
            % t = Nx1 or 1xN array
            assert(min(t)>=this.I(1),'t is out of I');
            assert(max(t)<=this.I(2),'t is out of I');
            maxDiff = max(diff([this.I(1);t(:)]));
            if maxDiff/diff(this.I)> 1/100
                filled = true;
                sizeSave = length(t);
                t = [t;linspace(this.I(1),this.I(2),1000)'];
                [t,indSort] = sort(t);
            end
            t_size = size(t);
            t = [this.I(1);t(:)];
            [xx,ww] = domGauss_Legendre1D(5,0,1);
            X = t(1:(end-1))*ones(1,length(xx)) + diff(t)*xx';
            W = diff(t)*ww';
            F = @(u)(sqrt(this.dx(u).^2 + this.dy(u).^2));
            Gs = F(X).*W;
            Gs = sum(Gs,2);
            Gs = cumsum(reshape(Gs,t_size(1),t_size(2)));
            if filled
                GS = Gs(1:sizeSave);
                for i = 1:sizeSave
                    GS(i) = Gs(indSort==i);
                end
                Gs = GS;
            end
        end
        function[t] = t_of_s(this,S)
            if this.isnormal
                t = S;
            else
                Sold = S;
                [S,p1] = sort(Sold);
                p2 = invPerm(p1);
                F = @(u)(sqrt(this.dx(u).^2 + this.dy(u).^2));
                size_s = size(S);
                S = S(:);
                tGuess = linspace(this.I(1),this.I(2),length(S));
                tGuess = tGuess(:);
                sGuess = s(this,tGuess);
                err = norm(sGuess - S,'inf');
                while err > 1e-6
                    tGuess = tGuess + (S - sGuess)./F(tGuess);
                    tGuess(tGuess < this.I(1)) = this.I(1);
                    tGuess(tGuess > this.I(2)) = this.I(2);
                    sGuess =  s(this,tGuess);
                    err = norm(sGuess - S,2);
                end
                t = tGuess;
                t = t(p2);
                t = reshape(t(:),size_s(1),size_s(2));
            end
        end
        function[c] = normalParam(this)
            X = this.x; Y = this.y;
            L = this.length;
            xx = @(s)(X(this.t_of_s(s)));
            yy = @(s)(Y(this.t_of_s(s)));
            ddx = this.dx;
            ddy = this.dy;
            c = SimpleCurve(xx,yy,[0,L],this.boundedSide);
            c = supplyDer(c,@(s)(auxDx(s)),@(s)(auxDy(s)));
            c.isnormal = true;
            function[dx] = auxDx(s)
                t = this.t_of_s(s);
                Dx = ddx(t); Dy = ddy(t);
                l = sqrt(Dx.^2 + Dy.^2);
                dx = Dx./l;
            end
            function[dy] = auxDy(s)
                t = this.t_of_s(s);
                Dx = ddx(t); Dy = ddy(t);
                l = sqrt(Dx.^2 + Dy.^2);
                dy = Dy./l;
            end
        end
        
        %% Curvature
        
        function[C] = curvature(this,t)
            t_size = size(t);
            t = t(:);
            dTx = safeDx(@Tx,this.I);
            dTy =  safeDx(@Ty,this.I);
            Txt = Tx(t);
            Tyt = Ty(t);
            dTxt = dTx(t(:));
            dTyt = dTy(t(:));
            num = Txt.*dTyt - Tyt.*dTxt;
            den = (Txt.^2 + Tyt.^2).^(3/2);
            C = num./den;
            if length(t) > 1                
                C = reshape(C,[t_size(1),t_size(2)]);
            end
            function[r] = Tx(t)
                [~,T] = this.tgt(t(:));
                r = T(:,1);
            end
            function[r] = Ty(t)
                [~,T] = this.tgt(t(:));
                r = T(:,2);
            end
            
        end
        
        
    end
end

