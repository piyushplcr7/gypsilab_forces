function [I] = specialIntegral2(alpha,a,b,X,N)

I = specialIntegral1(alpha,b,X,N) - specialIntegral2(alpha,a,X,N);

end



