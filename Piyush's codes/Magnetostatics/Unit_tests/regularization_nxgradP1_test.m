    clear; clc; close all;
    
    ivals = 5:5;
    Nivals = size(ivals,2);
    norms = zeros(Nivals,3);
    
    for i = 1:Nivals
        N = 2^ivals(i);
        mesh = mshSphere(N,1);
        kernel = cell(3,1);
        kernel{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
        kernel{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
        kernel{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);
        
        RWG = fem(mesh,'RWG');
        P1 = fem(mesh,'P1');
        DIV0 = nxgrad(P1);
        Gamma = dom(mesh,3);
        
        % mat 1 uses the regularization implemented by me
        mat1 = integral(Gamma,Gamma,RWG,kernel,DIV0);
        mat1 = mat1 + regularize(Gamma,Gamma,RWG,'grady[1/r]',DIV0);
        
        mat2 = integral(Gamma,Gamma,DIV0,kernel,RWG);
        mat2 = mat2 + regularize(Gamma,Gamma,DIV0,'grady[1/r]',RWG);

        norms(i,1) = norm(mat1);
        norms(i,2) = norm(mat2);
        norms(i,3) = norm(mat1-mat2');
    end

    plot(ivals,norms(:,1));
    hold on;
    plot(ivals,norms(:,2));
    plot(ivals,norms(:,3));
