
I         = gridd.I;
J         = gridd.J;
Nz        = gridd.Nz;
Ned       = gridd.Ned;
M         = gridd.M; % total number of nodes
MMe       = gridd.MMe;

%%%%% Pre-allocate some memory
Vbf       = zeros(I,Nz);
Vbb       = zeros(I,Nz);
c         = zeros(I,Nz);