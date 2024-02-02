function [Vel,DVel] = getTransVelDVelCutoff(dir,R)
    cutoff = @(X) vecnorm(X,2,2) < R;
    Vel = @(X) cutoff(X).*ones(size(X,1),1)* dir ;

    

    DVelmat = zeros(3,3);
    
    DVel = cell(3,1);

    DVel{1} = @(X) cutoff(X).*ones(size(X,1),1)* DVelmat(1,:);
    DVel{2} = @(X) cutoff(X).*ones(size(X,1),1)* DVelmat(2,:);
    DVel{3} = @(X) cutoff(X).*ones(size(X,1),1)* DVelmat(3,:);
end