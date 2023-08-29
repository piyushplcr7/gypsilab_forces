for i = 1:81
    disp(i);
    [shape_derivatives_bem(:,i) shape_derivatives_mst(:,i)]


end

% List of i values for which we see convergence
isphere = [1 4 7 10 13 16 19 22 25 28 29 30 37 38 39 46 47 48 55 56 57 58 59 60 61 62 63]

icube = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 37 38 39 46 47 48 55 56 57 60 61 62 ]

icuboid5 = [1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22 23 24 25 26 27 28 29 30 37 38 39 46 47 48 55 56 57 58 59 60 61 62 63]