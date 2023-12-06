function mshh = tetra3()
%  Matlab mesh
% untitled, Created by Gmsh
% ASCII
clear mshh;
mshh.nbNod = 130;
mshh.POS = [
1 0 -0.70710678118;
-1 0 -0.70710678118;
0 1 0.70710678118;
0 -1 0.70710678118;
-0.875000000000286 0.1249999999997141 -0.5303300858854044;
-0.7500000000006223 0.2499999999993777 -0.35355339059088;
-0.6250000000010584 0.3749999999989416 -0.1767766952964968;
-0.5000000000014828 0.4999999999985172 -2.096989248911996e-12;
-0.3750000000010666 0.6249999999989334 0.1767766952934916;
-0.2500000000006546 0.7499999999993454 0.3535533905890742;
-0.1250000000002899 0.8749999999997101 0.5303300858845901;
0.1249999999997141 0.875000000000286 0.5303300858854044;
0.2499999999993777 0.7500000000006223 0.35355339059088;
0.3749999999989416 0.6250000000010584 0.1767766952964968;
0.4999999999985172 0.5000000000014828 2.096989248911996e-12;
0.6249999999989334 0.3750000000010666 -0.1767766952934916;
0.7499999999993454 0.2500000000006546 -0.3535533905890742;
0.8749999999997101 0.1250000000002899 -0.5303300858845901;
0.875000000000286 -0.1249999999997141 -0.5303300858854044;
0.7500000000006223 -0.2499999999993777 -0.35355339059088;
0.6250000000010584 -0.3749999999989416 -0.1767766952964968;
0.5000000000014828 -0.4999999999985172 -2.096989248911996e-12;
0.3750000000010666 -0.6249999999989334 0.1767766952934916;
0.2500000000006546 -0.7499999999993454 0.3535533905890742;
0.1250000000002899 -0.8749999999997101 0.5303300858845901;
-0.1249999999997141 -0.875000000000286 0.5303300858854044;
-0.2499999999993777 -0.7500000000006223 0.35355339059088;
-0.3749999999989416 -0.6250000000010584 0.1767766952964968;
-0.4999999999985172 -0.5000000000014828 2.096989248911996e-12;
-0.6249999999989334 -0.3750000000010666 -0.1767766952934916;
-0.7499999999993454 -0.2500000000006546 -0.3535533905890742;
-0.8749999999997101 -0.1250000000002899 -0.5303300858845901;
-0.7500000000004167 0 -0.70710678118;
-0.5000000000011573 0 -0.70710678118;
-0.250000000001955 0 -0.70710678118;
-2.752797989558076e-12 0 -0.70710678118;
0.2499999999979354 0 -0.70710678118;
0.4999999999986235 0 -0.70710678118;
0.7499999999993117 0 -0.70710678118;
0 -0.7500000000004167 0.70710678118;
0 -0.5000000000011573 0.70710678118;
0 -0.250000000001955 0.70710678118;
0 -2.752797989558076e-12 0.70710678118;
0 0.2499999999979354 0.70710678118;
0 0.4999999999986235 0.70710678118;
0 0.7499999999993117 0.70710678118;
-0.2500000000007414 0.2499999999978822 0.3535533905889515;
-0.5 -1.482758360538128e-12 -1.308410382229149e-17;
-0.2499999999992586 -0.2500000000021178 0.3535533905910485;
-0.1250000000003707 0.1249999999975647 0.5303300858844758;
-0.1250000000003707 0.3749999999982529 0.5303300858844758;
-0.250000000000698 0.4999999999986138 0.3535533905890129;
-0.1250000000003273 0.6249999999989845 0.5303300858845371;
-0.3750000000011121 0.3749999999981998 0.1767766952934273;
-0.5000000000007414 0.2499999999985172 -1.048505539655682e-12;
-0.6250000000003111 0.1249999999989475 -0.17677669529544;
-0.6249999999996727 -0.1250000000010687 -0.1767766952945372;
-0.7499999999999838 -6.384337503106963e-13 -0.3535533905899772;
-0.4999999999992586 -0.2500000000014828 1.048479371448037e-12;
-0.3750000000003707 0.1249999999981997 0.1767766952944758;
-0.25 -2.117750419472486e-12 0.35355339059;
-0.3749999999996293 -0.1250000000018003 0.1767766952955243;
-0.1249999999996293 -0.1250000000024353 0.5303300858855243;
-0.3749999999988879 -0.3750000000018003 0.1767766952965728;
-0.2499999999993182 -0.50000000000137 0.3535533905909643;
-0.1249999999996293 -0.3750000000016376 0.5303300858855243;
-0.1249999999996888 -0.6250000000008897 0.5303300858854401;
0.5000000000000001 1.482758360538128e-12 3.925231146720425e-17;
0.2500000000007414 -0.250000000000635 0.3535533905889516;
0.2499999999992586 0.249999999999365 0.3535533905910486;
0.5000000000007415 -0.2499999999985172 -1.048453203240392e-12;
0.6250000000003113 -0.1249999999989475 -0.17677669529544;
0.6249999999996728 0.1250000000010687 -0.1767766952945371;
0.749999999999984 6.384337503106963e-13 -0.3535533905899771;
0.4999999999992587 0.2500000000014828 1.048531707863327e-12;
0.3749999999996294 0.1250000000004239 0.1767766952955243;
0.3749999999988879 0.3750000000004239 0.1767766952965728;
0.3750000000003707 -0.1249999999995761 0.1767766952944758;
0.25 -6.350475700855895e-13 0.3535533905900001;
0.3750000000011121 -0.3749999999995761 0.1767766952934273;
0.1250000000003707 -0.1250000000016939 0.5303300858844758;
0.1249999999996292 0.1249999999983061 0.5303300858855244;
0.2499999999993182 0.4999999999999936 0.3535533905909644;
0.1249999999996292 0.3749999999989942 0.5303300858855244;
0.1249999999996889 0.624999999999623 0.5303300858854401;
0.1250000000003707 -0.3750000000008962 0.5303300858844758;
0.250000000000698 -0.4999999999999902 0.353553390589013;
0.1250000000003273 -0.6250000000002514 0.5303300858845372;
0.2499999999978822 0.2500000000007414 -0.3535533905889516;
-1.482813871689359e-12 0.5000000000000001 -3.925231146720425e-17;
-0.2500000000021178 0.2499999999992586 -0.3535533905910486;
0.3749999999981997 0.3750000000011121 -0.1767766952934273;
0.4999999999986138 0.250000000000698 -0.353553390589013;
0.3749999999982528 0.1250000000003707 -0.5303300858844758;
0.6249999999989845 0.1250000000003273 -0.5303300858845372;
0.1249999999975647 0.1250000000003707 -0.5303300858844758;
-0.2500000000014828 0.4999999999992587 -1.048531707863327e-12;
-0.1250000000010687 0.6249999999996728 0.1767766952945371;
0.1249999999989475 0.6250000000003113 0.17677669529544;
-6.384337503106963e-13 0.749999999999984 0.3535533905899771;
0.2499999999985172 0.5000000000007415 1.048453203240392e-12;
0.1249999999981997 0.3750000000003707 -0.1767766952944758;
-0.1250000000018003 0.3749999999996294 -0.1767766952955243;
-2.117805930623717e-12 0.25 -0.3535533905900001;
-0.3750000000018003 0.3749999999988879 -0.1767766952965728;
-0.1250000000024353 0.1249999999996292 -0.5303300858855244;
-0.3750000000016376 0.1249999999996292 -0.5303300858855244;
-0.50000000000137 0.2499999999993182 -0.3535533905909644;
-0.6250000000008897 0.1249999999996889 -0.5303300858854401;
1.482758360538128e-12 -0.5000000000000001 -3.925231146753402e-17;
0.249999999999365 -0.2499999999992587 -0.3535533905910485;
-0.250000000000635 -0.2500000000007415 -0.3535533905889515;
-0.2499999999985172 -0.5000000000007415 1.048474061618001e-12;
-0.1249999999989475 -0.6250000000003112 0.1767766952954399;
0.1250000000010687 -0.6249999999996727 0.176776695294537;
6.384450141611725e-13 -0.7499999999999839 0.353553390589977;
0.2500000000014828 -0.4999999999992586 -1.048552566240936e-12;
0.1249999999983061 -0.1249999999996294 -0.5303300858855242;
-6.349755672092188e-13 -0.2500000000000001 -0.35355339059;
-0.1250000000016939 -0.1250000000003708 -0.5303300858844757;
0.1250000000004239 -0.3749999999996294 -0.1767766952955243;
-0.1249999999995761 -0.3750000000003708 -0.1767766952944758;
0.3750000000004239 -0.3749999999988881 -0.1767766952965728;
-0.3749999999995761 -0.3750000000011121 -0.1767766952934272;
0.4999999999999937 -0.2499999999993183 -0.3535533905909643;
0.3749999999989942 -0.1249999999996295 -0.5303300858855242;
0.624999999999623 -0.124999999999689 -0.5303300858854401;
-0.3750000000008962 -0.1250000000003707 -0.5303300858844757;
-0.4999999999999902 -0.250000000000698 -0.3535533905890128;
-0.6250000000002512 -0.1250000000003274 -0.530330085884537;
];
mshh.MAX = max(mshh.POS);
mshh.MIN = min(mshh.POS);
mshh.LINES =[
 2 5 0
 5 6 0
 6 7 0
 7 8 0
 8 9 0
 9 10 0
 10 11 0
 11 3 0
 3 12 0
 12 13 0
 13 14 0
 14 15 0
 15 16 0
 16 17 0
 17 18 0
 18 1 0
 1 19 0
 19 20 0
 20 21 0
 21 22 0
 22 23 0
 23 24 0
 24 25 0
 25 4 0
 4 26 0
 26 27 0
 27 28 0
 28 29 0
 29 30 0
 30 31 0
 31 32 0
 32 2 0
 2 33 0
 33 34 0
 34 35 0
 35 36 0
 36 37 0
 37 38 0
 38 39 0
 39 1 0
 4 40 0
 40 41 0
 41 42 0
 42 43 0
 43 44 0
 44 45 0
 45 46 0
 46 3 0
];
mshh.TRIANGLES =[
 43 50 44 0
 50 51 44 0
 50 47 51 0
 44 51 45 0
 47 52 51 0
 52 53 51 0
 52 10 53 0
 51 53 45 0
 47 54 52 0
 54 9 52 0
 54 8 9 0
 52 9 10 0
 45 53 46 0
 53 11 46 0
 53 10 11 0
 46 11 3 0
 8 55 7 0
 55 56 7 0
 55 48 56 0
 7 56 6 0
 48 57 56 0
 57 58 56 0
 57 31 58 0
 56 58 6 0
 48 59 57 0
 59 30 57 0
 59 29 30 0
 57 30 31 0
 6 58 5 0
 58 32 5 0
 58 31 32 0
 5 32 2 0
 8 54 55 0
 54 60 55 0
 54 47 60 0
 55 60 48 0
 47 61 60 0
 61 62 60 0
 61 49 62 0
 60 62 48 0
 47 50 61 0
 50 63 61 0
 50 43 63 0
 61 63 49 0
 48 62 59 0
 62 64 59 0
 62 49 64 0
 59 64 29 0
 29 64 28 0
 64 65 28 0
 64 49 65 0
 28 65 27 0
 49 66 65 0
 66 67 65 0
 66 41 67 0
 65 67 27 0
 49 63 66 0
 63 42 66 0
 63 43 42 0
 66 42 41 0
 27 67 26 0
 67 40 26 0
 67 41 40 0
 26 40 4 0
 22 71 21 0
 71 72 21 0
 71 68 72 0
 21 72 20 0
 68 73 72 0
 73 74 72 0
 73 17 74 0
 72 74 20 0
 68 75 73 0
 75 16 73 0
 75 15 16 0
 73 16 17 0
 20 74 19 0
 74 18 19 0
 74 17 18 0
 19 18 1 0
 15 75 77 0
 75 76 77 0
 75 68 76 0
 77 76 70 0
 68 78 76 0
 78 79 76 0
 78 69 79 0
 76 79 70 0
 68 71 78 0
 71 80 78 0
 71 22 80 0
 78 80 69 0
 70 79 82 0
 79 81 82 0
 79 69 81 0
 82 81 43 0
 15 77 14 0
 77 83 14 0
 77 70 83 0
 14 83 13 0
 70 84 83 0
 84 85 83 0
 84 45 85 0
 83 85 13 0
 70 82 84 0
 82 44 84 0
 82 43 44 0
 84 44 45 0
 13 85 12 0
 85 46 12 0
 85 45 46 0
 12 46 3 0
 43 81 42 0
 81 86 42 0
 81 69 86 0
 42 86 41 0
 69 87 86 0
 87 88 86 0
 87 24 88 0
 86 88 41 0
 69 80 87 0
 80 23 87 0
 80 22 23 0
 87 23 24 0
 41 88 40 0
 88 25 40 0
 88 24 25 0
 40 25 4 0
 15 92 16 0
 92 93 16 0
 92 89 93 0
 16 93 17 0
 89 94 93 0
 94 95 93 0
 94 38 95 0
 93 95 17 0
 89 96 94 0
 96 37 94 0
 96 36 37 0
 94 37 38 0
 17 95 18 0
 95 39 18 0
 95 38 39 0
 18 39 1 0
 8 97 9 0
 97 98 9 0
 97 90 98 0
 9 98 10 0
 90 99 98 0
 99 100 98 0
 99 13 100 0
 98 100 10 0
 90 101 99 0
 101 14 99 0
 101 15 14 0
 99 14 13 0
 10 100 11 0
 100 12 11 0
 100 13 12 0
 11 12 3 0
 15 101 92 0
 101 102 92 0
 101 90 102 0
 92 102 89 0
 90 103 102 0
 103 104 102 0
 103 91 104 0
 102 104 89 0
 90 97 103 0
 97 105 103 0
 97 8 105 0
 103 105 91 0
 89 104 96 0
 104 106 96 0
 104 91 106 0
 96 106 36 0
 36 106 35 0
 106 107 35 0
 106 91 107 0
 35 107 34 0
 91 108 107 0
 108 109 107 0
 108 6 109 0
 107 109 34 0
 91 105 108 0
 105 7 108 0
 105 8 7 0
 108 7 6 0
 34 109 33 0
 109 5 33 0
 109 6 5 0
 33 5 2 0
 29 113 28 0
 113 114 28 0
 113 110 114 0
 28 114 27 0
 110 115 114 0
 115 116 114 0
 115 24 116 0
 114 116 27 0
 110 117 115 0
 117 23 115 0
 117 22 23 0
 115 23 24 0
 27 116 26 0
 116 25 26 0
 116 24 25 0
 26 25 4 0
 36 118 120 0
 118 119 120 0
 118 111 119 0
 120 119 112 0
 111 121 119 0
 121 122 119 0
 121 110 122 0
 119 122 112 0
 111 123 121 0
 123 117 121 0
 123 22 117 0
 121 117 110 0
 112 122 124 0
 122 113 124 0
 122 110 113 0
 124 113 29 0
 22 123 21 0
 123 125 21 0
 123 111 125 0
 21 125 20 0
 111 126 125 0
 126 127 125 0
 126 38 127 0
 125 127 20 0
 111 118 126 0
 118 37 126 0
 118 36 37 0
 126 37 38 0
 20 127 19 0
 127 39 19 0
 127 38 39 0
 19 39 1 0
 36 120 35 0
 120 128 35 0
 120 112 128 0
 35 128 34 0
 112 129 128 0
 129 130 128 0
 129 31 130 0
 128 130 34 0
 112 124 129 0
 124 30 129 0
 124 29 30 0
 129 30 31 0
 34 130 33 0
 130 32 33 0
 130 31 32 0
 33 32 2 0
];
mshh.PNT =[
 1 0
 2 0
 3 0
 4 0
];

end
