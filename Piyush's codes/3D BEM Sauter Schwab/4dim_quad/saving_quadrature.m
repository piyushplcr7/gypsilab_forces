X = [];
W = [];
counter = 0;
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                counter = counter+1;
                X = [X; x(i) x(j) x(k) x(l)];
                W = [W; w(i) * w(j) * w(k) * w(l)];
            end
        end
    end
end