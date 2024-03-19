% Script to process data and compute dualnorm

function [conv_rate_mst_VP,conv_rate_bem_VP] = combineDualNormData(fname_VP,fname_SP,idxmst,idxbem,trim)
    load(fname_VP);

    shape_derivatives_bem_VP = shape_derivatives_bem(trim,:);
    shape_derivatives_mst_VP = shape_derivatives_mst(trim,:);
    hvals_VP = hvals(trim);

    N_VP = size(shape_derivatives_bem_VP,1);

    % Converged values chosen as the last value from shape_derivatives_bem
    converged_values = shape_derivatives_bem_VP(N_VP,:);

    errs_mst_VP = shape_derivatives_mst_VP - converged_values;
    errs_bem_VP = shape_derivatives_bem_VP - converged_values;

    % Conversion to dual norm
    gramat = gramat_abc_alpha();
    diagonal = diag(gramat)';
    diagonal = 1./diagonal;

    Ginv_errs_mst_VP = diagonal.*errs_mst_VP;
    Ginv_errs_bem_VP = diagonal.*errs_bem_VP;

    dualnorm_errs_mst_VP = sqrt(dot(errs_mst_VP,Ginv_errs_mst_VP,2));
    dualnorm_errs_bem_VP = sqrt(dot(errs_bem_VP,Ginv_errs_bem_VP,2));

    % Dualnorm errors SP
    load(fname_SP);

    shape_derivatives_bem_SP = shape_derivatives_bem(trim,:);
    shape_derivatives_mst_SP = shape_derivatives_mst(trim,:);
    hvals_SP = hvals(trim);

    N_SP = size(shape_derivatives_bem_SP,1);

    errs_mst_SP = shape_derivatives_mst_SP - converged_values;
    errs_bem_SP = shape_derivatives_bem_SP - converged_values;

    Ginv_errs_mst_SP = diagonal.*errs_mst_SP;
    Ginv_errs_bem_SP = diagonal.*errs_bem_SP;

    dualnorm_errs_mst_SP = sqrt(dot(errs_mst_SP,Ginv_errs_mst_SP,2));
    dualnorm_errs_bem_SP = sqrt(dot(errs_bem_SP,Ginv_errs_bem_SP,2));

    % Plotting

    figure;

    loglog(hvals_VP,dualnorm_errs_bem_VP,'-*', 'LineWidth', 2, 'MarkerSize', 10);
    hold on;
    loglog(hvals_VP,dualnorm_errs_mst_VP,'-o', 'LineWidth', 2, 'MarkerSize', 10);

    loglog(hvals_SP,dualnorm_errs_bem_SP,'-^', 'LineWidth', 2, 'MarkerSize', 10);
    loglog(hvals_SP,dualnorm_errs_mst_SP,'->', 'LineWidth', 2, 'MarkerSize', 10);
    
    hvals_small_VP = hvals_VP(1:N_VP-1);

    model_mst_VP = fitlm(log(hvals_VP(idxmst)),log(dualnorm_errs_mst_VP(idxmst)));
    model_bem_VP = fitlm(log(hvals_small_VP(idxbem)),log(dualnorm_errs_bem_VP(idxbem)));

    mst_coeffs_VP = model_mst_VP.Coefficients.Estimate;
    bem_coeffs_VP = model_bem_VP.Coefficients.Estimate;

    err_mst_fitted_VP = hvals_VP(idxmst).^mst_coeffs_VP(2) * exp(mst_coeffs_VP(1));
    err_bem_fitted_VP = hvals_small_VP(idxbem).^bem_coeffs_VP(2) * exp(bem_coeffs_VP(1));

    loglog(hvals_VP(idxmst),err_mst_fitted_VP,'--','Color',[0.7 0.7 0.7]);
    loglog(hvals_small_VP(idxbem),err_bem_fitted_VP,'--','Color',[0.7 0.7 0.7]);

    hvals_small_SP = hvals_SP(1:N_SP-1);

    model_mst_SP = fitlm(log(hvals_SP(idxmst)),log(dualnorm_errs_mst_SP(idxmst)));
    model_bem_SP = fitlm(log(hvals_small_SP(idxbem)),log(dualnorm_errs_bem_SP(idxbem)));

    mst_coeffs_SP = model_mst_SP.Coefficients.Estimate;
    bem_coeffs_SP = model_bem_SP.Coefficients.Estimate;

    err_mst_fitted_SP = hvals_SP(idxmst).^mst_coeffs_SP(2) * exp(mst_coeffs_SP(1));
    err_bem_fitted_SP = hvals_small_SP(idxbem).^bem_coeffs_SP(2) * exp(bem_coeffs_SP(1));

    loglog(hvals_SP(idxmst),err_mst_fitted_SP,'--','Color',[0.7 0.7 0.7]);
    loglog(hvals_small_SP(idxbem),err_bem_fitted_SP,'--','Color',[0.7 0.7 0.7]);
    
    % title(fname, 'Interpreter', 'none');

    xlabel('meshwidth');
    ylabel('dual norm error');

    conv_rate_mst_VP = mst_coeffs_VP(2);
    conv_rate_bem_VP = bem_coeffs_VP(2);
    conv_rate_mst_SP = mst_coeffs_SP(2);
    conv_rate_bem_SP = bem_coeffs_SP(2);

    legend(['BEM VP: ',num2str(conv_rate_bem_VP)],['MST VP: ', num2str(conv_rate_mst_VP)],...
        ['BEM SP: ',num2str(conv_rate_bem_SP)],['MST SP: ', num2str(conv_rate_mst_SP)],'Location','southeast');
    print([fname_VP(1:end-3),'eps'],'-depsc2');

end

