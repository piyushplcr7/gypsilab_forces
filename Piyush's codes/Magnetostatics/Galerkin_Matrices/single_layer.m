%  Single Layer BIO magnetostatics

function MV = single_layer(Gamma,test,trial)
        Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); 
        MV = integral(Gamma,Gamma,test,Gxy,trial)/(4*pi);
        MV = MV + 1/(4*pi)*regularize(Gamma,Gamma,test,'[1/r]',trial);
end

