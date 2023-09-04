function [Vel,DVel] = getTransVelDVel(dir)
    Vel = @(X) ones(size(X,1),1)* dir ;

    DVelmat = zeros(3,3);
    
    DVel = cell(3,1);

    DVel{1} = @(X) ones(size(X,1),1)* DVelmat(1,:);
    DVel{2} = @(X) ones(size(X,1),1)* DVelmat(2,:);
    DVel{3} = @(X) ones(size(X,1),1)* DVelmat(3,:);
end