% Double layer

function MK = double_layer(Gamma,test,trial)
    kernel = cell(3,1);
    kernel{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
    kernel{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
    kernel{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);

    % One more sign change due to flipped trial and test
    MK = 1/4/pi*integral(Gamma,Gamma,trial.nx,kernel,test);
    MK = MK + 1/4/pi*regularize(Gamma,Gamma,trial.nx,'grady[1/r]',test);

    % Transposing the matrix
    MK = -MK';
end