function [W,XX,YY] = getNestedQuadrature(test_gamma,trial_gamma)
    
    [X_test,W_test] = test_gamma.qud;
    [Y_trial,W_trial] = trial_gamma.qud;
    
    assert(size(X_test,2)==3);
    assert(size(Y_trial,2)==3);

    % repmat for trial, repelem for test
    NX = size(X_test,1);
    NY = size(T_trial,1);

    W = repelem(W_test,NY,1).*repmat(W_trial,NX,1);
    XX = repelem(X_test,NY,1);
    YY = repmat(Y_trial,NX,1);

end