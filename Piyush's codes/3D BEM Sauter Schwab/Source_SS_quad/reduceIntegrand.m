function integral = reduceIntegrand(integrand,NI,Nqud)
    NI0 = NI(1);
    NI1 = NI(2);
    NI2 = NI(3);
    NI3 = NI(4);
    NQ0 = Nqud(1);
    NQ1 = Nqud(2);
    NQ2 = Nqud(3);
    NQ3 = Nqud(4);

    dim = size(integrand,2);
    integrand1 = integrand(1:NI0*NQ0,:);
    integrand2 = integrand(NI0*NQ0+1:NI0*NQ0+NI1*NQ1,:);
    integrand3 = integrand(NI0*NQ0+NI1*NQ1+1:NI0*NQ0+NI1*NQ1+NI2*NQ2,:);
    integrand4 = integrand(NI0*NQ0+NI1*NQ1+NI2*NQ2+1:end,:);


    integral = [reshape(sum(reshape(integrand1,NQ0,NI0*dim),1),NI0,dim);...
                reshape(sum(reshape(integrand2,NQ1,NI1*dim),1),NI1,dim);...
                reshape(sum(reshape(integrand3,NQ2,NI2*dim),1),NI2,dim);...
                reshape(sum(reshape(integrand4,NQ3,NI3*dim),1),NI3,dim)];
end