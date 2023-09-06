% Gramian matrix

addpath(genpath("../../../"));
clear; clc; close all;
%format long;

% side
l = 6*pi;

N = 5000;

mesh = mshCube(N,[l l l]);

omega = dom(mesh,4);

n = 3;

manually_constructed = zeros(n^3,1);

% Compute the gramian matrix. No. of basis functions = 3*n^3
gramatfull = zeros(3*n^3,3*n^3);
gramat_small = zeros(n^3,n^3);

% With zero indexing, basis[i+j*n+k*n^2] = getSinVelDVel(i,j,k);

for II = 1:n^3
    % convert II,JJ to a basis
    i1 = mod(II-1,n)+1;
    k1 = floor((II-1)/n^2)+1;
    j1 = floor((II-1)/n);
    j1 = mod(j1,n)+1;
    %j1 = j1-(k1-1)*n^2+1;
    manually_constructed(II) = (l/2)^3 * (1 + i1^2 + j1^2 + k1^2);
    %disp(II);
    fprintf("%d Corresponds to basis %d %d %d \n",II,i1,j1,k1);
    for JJ = 1:n^3
        

        i2 = mod(JJ-1,n)+1;
        k2 = floor((JJ-1)/n^2)+1;
        j2 = floor((JJ-1)/n);
        j2 = mod(j2,n)+1;
        %j2 = j2-(k2-1)*n^2+1;

        [Vel1,DVel1] = getCosVelDVel(i1,j1,k1,1);
        [Vel2,DVel2] = getCosVelDVel(i2,j2,k2,1);

%         [Vel1,DVel1] = getSinVelDVel(i1,j1,k1,1);
%         [Vel2,DVel2] = getSinVelDVel(i2,j2,k2,1);

        [X,W] = omega.qud;

        % L2 inner pdt
        l2pdt = sum(W.*dot(Vel1(X),Vel2(X),2),1);

        % H1 inner pdt
        h1pdt = sum(W.*(dot(DVel1{1}(X),DVel2{1}(X),2)...
                       +dot(DVel1{2}(X),DVel2{2}(X),2)...
                       +dot(DVel1{3}(X),DVel2{3}(X),2)),1);

        % h1pdt = sum(W.*dot(DVel1{1}(X),DVel2{1}(X),2),1);

        gramat_small(II,JJ) = l2pdt + h1pdt;

    end
end

gramatfull(1:n^3,1:n^3) = gramat_small;
gramatfull(n^3+1:2*n^3,n^3+1:2*n^3) = gramat_small;
gramatfull(2*n^3+1:3*n^3,2*n^3+1:3*n^3) = gramat_small;

% for dir = 1:3
%     for i = 1:n
%         for j = 1:n
%             for k = 1:n
%                 [Vel,DVel] = getSinVelDVel(i,j,k,dir);
%             end
%         end
%     end
% end



