
I          = gridd.I ;
%%%%% Pre-allocate some memory
dVf        = zeros(I,Nz);
dVb        = zeros(I,Nz);
c          = zeros(I,Nz);
s          = param.s;         % CRRA
rho        = param.rho;       % Subjective discount rate
r          = guess.r;         % Interest rate
la1        = param.la1;       % Transition rates
la2        = param.la2;
la3        = param.la3;       % Transition rates
la4        = param.la4;
la5        = param.la5;
la6        = param.la6;
lx1        = param.lx1;
lx2        = param.lx2;
lx3        = param.lx3;
constant   = param.constant;