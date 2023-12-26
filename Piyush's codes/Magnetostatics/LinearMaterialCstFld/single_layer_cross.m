%  Single Layer BIO magnetostatics cross
%  Kernel should not hit the singularity!

function MV = single_layer_cross(Gamma_test, Gamma_trial,test,trial)
        Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); 
        MV = integral(Gamma_test,Gamma_trial,test,Gxy,trial)/(4*pi);
end

