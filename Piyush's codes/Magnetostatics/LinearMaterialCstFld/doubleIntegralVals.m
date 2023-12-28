function integral = doubleIntegralVals(Gamma_test,Gamma_trial,test_vals,kernel,trial_vals)
    [X,WX] = Gamma_test.qud;
    [Y,WY] = Gamma_trial.qud;
    NX = size(X,1);
    NY = size(Y,1);
    XX = repmat(X,NY,1); WWX = repmat(WX,NY,1); 
    YY = repelem(Y,NX,1); WWY = repelem(WY,NX,1);
    W = WWX .* WWY;

    integral =  sum(W.*dot)

end