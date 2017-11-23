function [A,B] = hmxLowrank(Mh)
%+========================================================================+
%|                                                                        |
%|         OPENHMX - LIBRARY FOR H-MATRIX COMPRESSION AND ALGEBRA         |
%|           openHmx is part of the GYPSILAB toolbox for Matlab           |
%|                                                                        |
%| COPYRIGHT : Matthieu Aussal (c) 2015-2017.                             |
%| PROPERTY  : Centre de Mathematiques Appliquees, Ecole polytechnique,   |
%| route de Saclay, 91128 Palaiseau, France. All rights reserved.         |
%| LICENCE   : This program is free software, distributed in the hope that|
%| it will be useful, but WITHOUT ANY WARRANTY. Natively, you can use,    |
%| redistribute and/or modify it under the terms of the GNU General Public|
%| License, as published by the Free Software Foundation (version 3 or    |
%| later,  http://www.gnu.org/licenses). For private use, dual licencing  |
%| is available, please contact us to activate a "pay for remove" option. |
%| CONTACT   : matthieu.aussal@polytechnique.edu                          |
%| WEBSITE   : www.cmap.polytechnique.fr/~aussal/gypsilab                 |
%|                                                                        |
%| Please acknowledge the gypsilab toolbox in programs or publications in |
%| which you use it.                                                      |
%|________________________________________________________________________|
%|   '&`   |                                                              |
%|    #    |   FILE       : hmxLowrank.m                                  |
%|    #    |   VERSION    : 0.31                                          |
%|   _#_   |   AUTHOR(S)  : Matthieu Aussal                               |
%|  ( # )  |   CREATION   : 14.03.2017                                    |
%|  / 0 \  |   LAST MODIF : 25.11.2017                                    |
%| ( === ) |   SYNOPSIS   : Convert H-Matrix to low-rank approximation    |
%|  `---'  |                                                              |
%+========================================================================+

% H-Matrix (recursion)
if (Mh.typ == 0)
    % Initialisation
    A = zeros(Mh.dim(1),0);
    B = zeros(0,Mh.dim(2));

    % Recursion
    for i = 1:4
        % Children evaluation
        [Ai,Bi] = hmxLowrank(Mh.chd{i});

        % Low-rank addition for A
        tmp              = zeros(Mh.dim(1),size(Ai,2));
        tmp(Mh.row{i},:) = Ai;
        A                = [A , tmp];

        % Low-rank addition for B
        tmp              = zeros(size(Bi,1),Mh.dim(2));
        tmp(:,Mh.col{i}) = Bi;
        B                = [B ; tmp];
    end

    % Recompression
    [A,B] = hmxQRSVD(A,B,Mh.tol);

% Compressed leaf
elseif (Mh.typ == 1)
    A = Mh.dat{1};
    B = Mh.dat{2};

% Full leaf
elseif (Mh.typ == 2)
    [A,B,flag] = hmxACA(Mh.dat,Mh.tol);
    if ~flag
        A = Mh.dat;
        B = eye(Mh.dim(2));
    end

% Others    
else
    error('hmxLowrank.m : unavailable case')
end
end
