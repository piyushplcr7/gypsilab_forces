% Checking something

% bndmesh = mshSphere(20,1);
bndmesh = mshCube(20,[1 1 1]);
bndmesh = bndmesh.bnd;

Gamma = dom(bndmesh,3);

NED = fem(bndmesh,'NED');
RWG = fem(bndmesh,'RWG');

Mnn = mass_matrix(Gamma,NED,NED);
Mrr = mass_matrix(Gamma,RWG,RWG);

Mnn_curl_curl = mass_matrix(Gamma,NED.curl,NED.curl);