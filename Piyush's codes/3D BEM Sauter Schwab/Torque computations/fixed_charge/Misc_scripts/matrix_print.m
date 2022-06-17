% Script to print the matrices in right format
% fprintf("[");
% for i = 1:81
%     fprintf("%.17d,",W(i));
% end
% fprintf("\n");

fprintf("[");
for i = 1:81
    for j = 1:4
        if j  == 4
            fprintf("%.17d;",X(i,j));
        else 
            fprintf("%.17d,",X(i,j));
        end
    end
    
end
fprintf("]\n");