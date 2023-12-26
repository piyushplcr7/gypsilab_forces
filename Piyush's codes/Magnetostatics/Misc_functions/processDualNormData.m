% Script to process data and compute dualnorm

function [conv_rate_mst,conv_rate_bem] = processDualNormData(fname,idxmst,idxbem,trim)
    load(fname);

    shape_derivatives_bem = shape_derivatives_bem(trim,:);
    shape_derivatives_mst = shape_derivatives_mst(trim,:);
    hvals = hvals(trim);

    N = size(shape_derivatives_bem,1);

    % Converged values chosen as the last value from shape_derivatives_bem
    converged_values = shape_derivatives_bem(N,:);

    errs_mst = shape_derivatives_mst - converged_values;
    errs_bem = shape_derivatives_bem(1:N-1,:) - converged_values;

    % Conversion to dual norm
    gramat = gramat_abc_alpha();
    diagonal = diag(gramat)';
    diagonal = 1./diagonal;

    Ginv_errs_mst = diagonal.*errs_mst;
    Ginv_errs_bem = diagonal.*errs_bem;

    dualnorm_errs_mst = sqrt(dot(errs_mst,Ginv_errs_mst,2));
    dualnorm_errs_bem = sqrt(dot(errs_bem,Ginv_errs_bem,2));

    figure;

    loglog(hvals,dualnorm_errs_mst,'-*');
    hold on;
    loglog(hvals(1:N-1),dualnorm_errs_bem,'-+');
    xlabel('h');
    ylabel('Dual norm error');
    
    hvals_small = hvals(1:N-1);

    model_mst = fitlm(log(hvals(idxmst)),log(dualnorm_errs_mst(idxmst)));
    model_bem = fitlm(log(hvals_small(idxbem)),log(dualnorm_errs_bem(idxbem)));

    mst_coeffs = model_mst.Coefficients.Estimate;
    bem_coeffs = model_bem.Coefficients.Estimate;

    err_mst_fitted = hvals(idxmst).^mst_coeffs(2) * exp(mst_coeffs(1));
    err_bem_fitted = hvals_small(idxbem).^bem_coeffs(2) * exp(bem_coeffs(1));

    loglog(hvals(idxmst),err_mst_fitted,'--','Color',[0.7 0.7 0.7]);
    loglog(hvals_small(idxbem),err_bem_fitted,'--','Color',[0.7 0.7 0.7]);
    
    title(fname(1:end-4), 'Interpreter', 'none');

    conv_rate_mst = mst_coeffs(2);
    conv_rate_bem = bem_coeffs(2);
    legend1 = ['VOL: ', num2str(conv_rate_mst)];
    l1 = legend1;
    legend2 = ['BEM: ', num2str(conv_rate_bem)];
    l2 = legend2;
    legend(l1,l2);
    print([fname(1:end-3), 'eps'], '-depsc2');

end

