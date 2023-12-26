% Double layer cross

function MK = double_layer_magnetostatics_cross(Gamma_test, Gamma_trial,test,trial)
    kernel = cell(3,1);
    kernel{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
    kernel{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
    kernel{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);

    % One more sign change due to flipped trial and test
    MK = 1/4/pi*integral(Gamma_test,Gamma_trial,test,kernel,trial);
%     MK = MK + 1/4/pi*regularize(Gamma_test,Gamma_test,test,'grady[1/r]',trial);
    MK = -MK;
end