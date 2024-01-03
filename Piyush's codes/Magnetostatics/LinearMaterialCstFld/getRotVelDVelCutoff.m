function [Vel,DVel] = getRotVelDVelCutoff(axis,Xcg,R)
    cutoff = @(X) vecnorm(X,2,2) < R;

    Vel = @(X) cutoff(X).*cross(ones(size(X,1),1)* axis ,X-Xcg);
    % Get the DVel matrix based on the axis a
    DVelmatt = @(a) [0 -a(3) a(2); a(3) 0 -a(1); -a(2) a(1) 0];
    DVelmat = DVelmatt(axis);
    
    DVel = cell(3,1);

    DVel{1} = @(X) cutoff(X).*ones(size(X,1),1)* DVelmat(1,:);
    DVel{2} = @(X) cutoff(X).*ones(size(X,1),1)* DVelmat(2,:);
    DVel{3} = @(X) cutoff(X).*ones(size(X,1),1)* DVelmat(3,:);
end