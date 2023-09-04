function DcurlA = compute_vecpot_D_curl(J,omega_src,X)
   % Integration happens over omega_src, variable called y
   [Ysrc,Wsrc] = omega_src.qud;

   Nsrc = size(Ysrc,1);
   Neval = size(X,1);

   X_ex = repelem(X,Nsrc,1);
   Ysrc_ex = repmat(Ysrc,Neval,1);

   xmy = X_ex-Ysrc_ex;

   % Rows of the matrix gradx gradx G for the evaluation points
   % X_ex,Ysrc_ex
   row1 = 3/(4*pi)* xmy(:,1).*xmy./vecnorm(xmy,2,2).^5 - 1/(4*pi) * [1 0 0]./vecnorm(xmy,2,2).^3; 
   row2 = 3/(4*pi)* xmy(:,2).*xmy./vecnorm(xmy,2,2).^5 - 1/(4*pi) * [0 1 0]./vecnorm(xmy,2,2).^3; 
   row3 = 3/(4*pi)* xmy(:,3).*xmy./vecnorm(xmy,2,2).^5 - 1/(4*pi) * [0 0 1]./vecnorm(xmy,2,2).^3; 

   % Evaluating the source current at the qudrature points
   J_Ysrc_ex = J(Ysrc_ex(:,1),Ysrc_ex(:,2),Ysrc_ex(:,3));

   % rows of the cross product
   row1 = cross(row1,J_Ysrc_ex,2);
   row2 = cross(row2,J_Ysrc_ex,2);
   row3 = cross(row3,J_Ysrc_ex,2);
   
   % Integration for row1
   Wsrc_ex = repmat(Wsrc,Neval,1);
   val1 = sum(reshape(Wsrc_ex .* row1,[Nsrc,Neval * 3]),1);

   % Integration for row2
   val2 = sum(reshape(Wsrc_ex .* row2,[Nsrc,Neval * 3]),1);

   % Integration for row3
   val3 = sum(reshape(Wsrc_ex .* row3,[Nsrc,Neval * 3]),1);

   % NEED TO TRANSPOSE TO GET D curl A because we compute grad curl A
   % rows -> cols
   row1 = reshape(val1,[Neval,3]); 
   row2 = reshape(val2,[Neval,3]);
   row3 = reshape(val3,[Neval,3]);

   DcurlA = cell(3,1);
   % Manually transposing
   DcurlA{1} = [row1(:,1) row2(:,1) row3(:,1)];
   DcurlA{2} = [row1(:,2) row2(:,2) row3(:,2)];
   DcurlA{3} = [row1(:,3) row2(:,3) row3(:,3)];


   % Other way to do these integrations
%    val1 = sum(Wsrc.*reshape(row1(:,1),[Nsrc,Neval]),1);
%    val2 = sum(Wsrc.*reshape(row1(:,2),[Nsrc,Neval]),1);
%    val3 = sum(Wsrc.*reshape(row1(:,3),[Nsrc,Neval]),1);

end