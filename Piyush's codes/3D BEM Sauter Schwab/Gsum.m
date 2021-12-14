function out = Gsum(n,alpha)

%     nidx = 0:n;
%     nidx2 = 2 * nidx;
%     Gvec = G(nidx', alpha * (nidx'==nidx'));
% 
%     out = sum(Gvec.*((2*nidx'+1).*Gvec - ) 

    out = 0;
    for i = 0:n
        out = out + G(i,alpha) * ( (2*i+1)*G(i,alpha) - (2*i+2)*G(i+1,alpha) );
    end
    
end