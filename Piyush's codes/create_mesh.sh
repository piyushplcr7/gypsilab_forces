#rm -rf results.txt

NN=20
g++ geom.cpp -I/usr/include/eigen3 -DMM=${NN}
./a.out

./gmsh -2 kite_sq.geo -o sqkite_mesh.msh -format msh2 -save_all
