function m = linfitandplot(h,err)
    mdl = fitlm(log(h),log(err));
    c = mdl.Coefficients.Estimate(1);
    m = mdl.Coefficients.Estimate(2);
    % Plotting the fitted line
    loglog(h,h.^m * exp(c),'--',Color="black");
end