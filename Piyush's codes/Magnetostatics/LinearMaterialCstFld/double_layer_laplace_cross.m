function MK = double_layer_laplace_cross(Gamma_test, Gamma_trial,test,trial)
    kernel = cell(3,1);
    kernel{1} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]1',0);
    kernel{2} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]2',0);
    kernel{3} = @(X,Y)femGreenKernel(X,Y,'grady[1/r]3',0);

    % One more sign change due to flipped trial and test
    MK = 1/4/pi*integral(Gamma_test,Gamma_trial,test,kernel,ntimes(trial));

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
%         kernel = 1/4/pi * (XX-YY)./vecnorm(XX-YY,2,2).^3;
% 
%         % Assuming vectorial test and trial spaces for this case
%         uqm_trial = trial.uqm(Gamma_trial);
% 
%         uqm_test = test.uqm(Gamma_test);
%         Ndof_trial = trial.ndof;
%         Ndof_test = test.ndof;
% 
%         normals_trial = Gamma_trial.qudNrm;
% 
%         mat = zeros(Ndof_test,Ndof_trial);
% 
%         for i = 1:Ndof_test
%             beta_i = uqm_test(:,i);
%             betaiXX = repelem(beta_i,N_trial,1);
%             for j = 1:Ndof_trial
%                 b_j = uqm_trial(:,j).*normals_trial;
%                 bjYY = repmat(b_j,N_test,1);
%                 mat(i,j) = sum(W.*dot(kernel,bjYY,2).*betaiXX,1);
%             end
% 
%         end
%         disp('double layer cross check');
%         norm(mat-MK)
end