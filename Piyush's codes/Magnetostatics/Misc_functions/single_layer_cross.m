%  Single Layer BIO magnetostatics cross
%  Kernel should not hit the singularity!

function MV = single_layer_cross(Gamma_test, Gamma_trial,test,trial)
        Gxy = @(X,Y)femGreenKernel(X,Y,'[1/r]',0); 
        MV = integral(Gamma_test,Gamma_trial,test,Gxy,trial)/(4*pi);

        % Alternative computation (explicitly)

%         [X_test,W_test] = Gamma_test.qud;
%         [Y_trial,W_trial] = Gamma_trial.qud;
%         N_test = size(X_test,1);
%         N_trial = size(Y_trial,1);
% 
%         XX = repelem(X_test,N_trial,1);
%         YY = repmat(Y_trial,N_test,1);
% 
%         W = repelem(W_test,N_trial,1).*repmat(W_trial,N_test,1);
% 
%         kernel = 1/4/pi./vecnorm(XX-YY,2,2);
% 
%         % Assuming vectorial test and trial spaces for this case
%         uqm_cell_trial = trial.uqm(Gamma_trial);
% 
%         uqm_cell_test = test.uqm(Gamma_test);
%         Ndof_trial = trial.ndof;
%         Ndof_test = test.ndof;
% 
%         mat = zeros(Ndof_test,Ndof_trial);
% 
%         for i = 1:Ndof_test
% %             beta_i = [uqm_cell_test{1}(:,i) uqm_cell_test{2}(:,i) uqm_cell_test{3}(:,i)];
%             beta_i = uqm_cell_test(:,i);
%             betaiXX = repelem(beta_i,N_trial,1);
%             for j = 1:Ndof_trial
% %                 b_j = [uqm_cell_trial{1}(:,j) uqm_cell_trial{2}(:,j) uqm_cell_trial{3}(:,j)];
%                 b_j = uqm_cell_trial(:,j);
%                 bjYY = repmat(b_j,N_test,1);
% %                 mat(i,j) = sum(W.*kernel.*dot(betaiXX,bjYY,2),1);
%                 mat(i,j) = sum(W.* kernel.* betaiXX.*bjYY,1);
%             end
% 
%         end
%         disp('single layer cross check');
%         norm(mat-MV)
end

