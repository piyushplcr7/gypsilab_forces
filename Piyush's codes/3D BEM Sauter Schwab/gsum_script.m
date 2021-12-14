% plot gsum
alpha = 1.05;
nvals = 1:100;
vals = nvals * 0;
for i = nvals
    vals(i) = Gsum(i,alpha);
end
plot(nvals,vals)