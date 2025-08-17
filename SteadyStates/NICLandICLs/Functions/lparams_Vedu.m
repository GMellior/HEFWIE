
crra       = param.s;         % CRRA
rho        = param.rho;       % Subjective discount rate
la1        = param.la1;       % Transition rates
la2        = param.la2;
la3        = param.la3;
la4        = param.la4;
constant   = param.constant;  % Constant in utility unction
ra         = param.ra;        % Interest and amortisation
r          = guess.r;
amort      = param.amort;

% Build b,a,income mesh
I         = gridd.I;
J         = gridd.J;
Nz        = gridd.Nz;
Ned       = gridd.Ned;
M         = gridd.M;          % Total number of nodes
MMe       = gridd.MMe;
amin      = gridd.amin;
amax      = gridd.amax;
bmin      = gridd.bmin;
bmax      = gridd.bmax;
a         = gridd.a;
aa        = gridd.aa;
bb        = gridd.bb;
aaa       = gridd.aaa;
bbb       = gridd.bbb;
da        = gridd.da;

% Numerical
tol_hjb   = num.tol_hjb;     % Convergence criteria
maxit     = num.maxit_hjb;   % Max number of loops
Delta     = num.Delta;       % Time step

%%%%% Pre-allocate some memory
Vbf       = zeros(I,J,Ned);
Vbb       = zeros(I,J,Ned);
c         = zeros(I,J,Ned);