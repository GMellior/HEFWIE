function Res = SSTSGT(subsidy,GT,Edpricescale,smallgrid)
% Solves for steady states (SS) of TS and GT
% Gustavo Mellior g.mellior@liverpool.ac.uk

% ***** INPUTS *****
% 1- Subsidy -  Edu cost % paid by government
% 2- TS or GT
% 3- Scale edu cost

warning('off','all')
Res.codename   = mfilename;                                                % Store filename
fniniguess     = 'startval';                                               % Load baseline parameters
load(fniniguess)                                                           % Use accross all steady states to keep things consistent with baseline
guess.r        =  0.05;                                                    % Guess interest rate
r              = guess.r;
xi             = xiBENCH;                                                  % CES share
z_ave          = z_aveBENCH;                                               % Effective labour supply
KD             = KBENCH;                                                   % Initial guess capital stock
Kold           = KD;
itax           = itaxBENCH;                                                % Initial guess of income tax rate
P              = PBENCH;                                                   % Initial guess of cost of education
param.altruism = altr;                                                     % Altruism strength
param.mcost    = mcBENCH*ones(1,2);                                        % Mental cost
mcost          = param.mcost;                                              
rho            = rhoBENCH;                                                 % Subjective discount rate
param.rho      = rho;
% Guess population shares (of workers) 
unEmp          = popBENCH(2);                                              % Worker without college degree
eEmp           = popBENCH(5);                                              % Worker with college degree
stn            = popBENCH(3);                                              % Student
cwp_enforced   = 0;                                                        % Enforce CWP 1.7 (only relevant for baseline, so keep at zero)
GDPscale       = (328.795/206.538)*54541;                                  % See Appendix
EY             = 0.79;
P_GDPscale     = (((0.75*21370+0.25*37430)*0.6 + 48510*0.4)/GDPscale)*Edpricescale;
crra           = 2;                                                        % Risk aversion in CRRA utility
param.s        = crra;                                                     % Store in structure used by functions
param.alpha    = 0.36;                                                     % Elasticity of Y wrt K
constant       = 0;                                                        % Constant in CRRA utility function
param.constant = constant;                                                 % Store in structure used by functions
constantkfe    = 1/4;                                                      % Inverse avg time until completion of university degree 1/constantkfe
constantkfe2   = 1000000;                                                  % For oldschool2.m
param.kappa    = 1/49;                                                     % Death rate
elastCES       = 2.5;                                                      % CES labour
nu             = (elastCES-1)/elastCES;
ipc_anddnm     = P_GDPscale*(1/constantkfe);                               % Edu cost per capita scale
Aprod          = 1;                                                        % Productivity constant
d              = dBENCH;                                                   % Depreciation
zstn           = 0.675*0.75;                                               % Labour efficiency of students
lx1            = 0.048;                                                    % Edu skills depreciation rate for unemployed
lx2            = 0.048/2;                                                  % Edu skills depreciation rate for employed
param.lx1      = lx1;                                                      %
param.lx2      = lx2;                                                      %
param.lx3_vec  = [1-((0.325+0.519)/2)^(1/4), 1-((0.207+0.256)/2)^(1/4)];   % College dropout risk
param.lx3      = param.lx3_vec;
parent_abtran  = [0.74136275, 0.25863725; 0.34289128, 0.65710872; 0.62180579, 0.37819421; 0.26403233, 0.7359677]; % Intergen trans of ability
param.la1      = 0.157*12;                                                 % Rate from unemp to emp   with no edu 
param.la2      = 0.011*12;                                                 % Rate from emp   to unemp with no edu 
param.la3      = 0.134*12;                                                 % Rate from unemp to emp   with edu 
param.la4      = 0.006*12;                                                 % Rate from emp   to unemp with edu
param.la5      = (1/3)*constantkfe;                                        % Graduating without a job
param.la6      = (2/3)*constantkfe;                                        % Graduating with    a job
reprate        = 0.382;                                                    % Unemployment benefits replacement rate
abtypes        = 2;                                                        % Ability types

%%

onerun                = tic;                                               % Track runtime
delta_kfe             = 1000;                                              % KFE update (unused)                                            
param.AssumeConcavity = 0;                                                 % Do not change - helps for non concavity of V
% Number of grid points
switch smallgrid
    case 1
        gridd.I       = 204;                                           
    otherwise
        if subsidy>0.7
            gridd.I   = 228;                                       
        else
            gridd.I   = 204;                                       
        end
end
g                     = linspace(0.024,0.039,gridd.I);                     % Growth rate of wealth step
new                   = (((1+g).^(1:gridd.I)'-1)./((1+g).^(gridd.I)'-1));
gridd.bmin            = -(20571/GDPscale)*EY;                              % Overdraft limit
gridd.bmax            = 50;                                                % Max wealth value
gridd.b               = [gridd.bmin; gridd.bmin + (gridd.bmax-gridd.bmin).*new];
gridd.I               = numel(gridd.b);                                    % Number of grid points stored in gridd
I                     = gridd.I;
gridd.bmin            = min(gridd.b(1),0);
gridd.bzeroind        = max(find(gridd.b<0));                              
% Income and income taxes
gridd.Nz              = 2;                                                 % Number of employment states without HE
gridd.Ned             = 3;                                                 % Number of employment states with HE + student type
Ned                   = gridd.Ned;
Nz                    = gridd.Nz;
NzNed                 = Nz+gridd.Ned;                                      % Total number of types
J                     = 1;                                                 % Grid points in "a" dimension (no student loans)
gridd.J               = J;
gridd.M               = I*J;                                               % Total number of nodes per type
M                     = gridd.M;
MMu                   = M*Nz;                                              % Total grid points without degree
gridd.MMu             = MMu;
gridd.MMe             = M*gridd.Ned;
MMe                   = gridd.MMe;                                         % Total gridd points with HE (including students)
MMM                   = MMu+MMe;                                           % Total grid points per ability type
MMMM                  = MMM*abtypes;                                       % Total grid points
% Intergenerational flows
parent_ab_stay        = param.kappa*[[parent_abtran(2,2)*ones(MMe,1);parent_abtran(4,2)*ones(MMu,1)],...
                                     [parent_abtran(1,1)*ones(MMe,1);parent_abtran(3,1)*ones(MMu,1)]]; 
parent_abs_small      = (1-[[parent_abtran(2,2);parent_abtran(4,2)],...
                        [parent_abtran(1,1);parent_abtran(3,1)]]);
parent_ab_stay2       = param.kappa*(1-[[parent_abtran(2,2)*ones(MMe,1);parent_abtran(4,2)*ones(MMu,1)],...
                                        [parent_abtran(1,1)*ones(MMe,1);parent_abtran(3,1)*ones(MMu,1)]]); 
% Initial guess for K demand, wages and taxes
gridd.wageL           = ((1-param.alpha)*Aprod)*(KD.^param.alpha)*(z_ave^(1-param.alpha-nu))*(xi*(unEmp+zstn*stn)^(nu-1));
gridd.wageH           = ((1-param.alpha)*Aprod)*(KD.^param.alpha)*(z_ave^(1-param.alpha-nu))*((1-xi)*(eEmp)^(nu-1));
wL                    = gridd.wageL;
wH                    = gridd.wageH;
gridd.income          = [gridd.wageL*reprate gridd.wageL];                 % Unemp and emp no college
gridd.income          = [gridd.income gridd.wageL*zstn gridd.wageH*reprate gridd.wageH]; 
gridd.iiincome        = zeros(gridd.I,NzNed);
param.itax2           = gridd.iiincome;
param.itax            = [0 itax 0 0 itax];                                 % Guess income tax
for i=1:NzNed
    gridd.iiincome(:,i)   = gridd.income(1,i);
    param.itax2(:,i)      = kron(param.itax(1,i),ones(gridd.I,J));
end
P0                    = (1-subsidy)*P;                                     % Education price individuals face

la                    = [param.la1,param.la2];
% Delta b (adjustments for nonuniform grids)
dbgrid                = diff(gridd.b);
gridd.dbgridf         = [dbgrid; dbgrid(I-1)];
gridd.dbgridb         = [dbgrid(1); dbgrid];
gridd.dbyF            = repmat(gridd.dbgridf(1:I-1),1,gridd.Ned); 
gridd.dbyB            = repmat(gridd.dbgridb(2:I),1,gridd.Ned);
% Trapezoidal rule: for KFE and moments
gridd.bdelta          = zeros(I,1);
gridd.bdelta(1)       = 0.5*dbgrid(1);
gridd.bdelta(2:I-1)   = 0.5*dbgrid(1:I-2) + 0.5*dbgrid(2:I-1);
gridd.bdelta(I)       = 0.5*dbgrid(I-1);
gridd.bdeltaNed       = repmat(gridd.bdelta,1,gridd.Ned);
baydelta              = repmat(gridd.bdelta,1,NzNed);                      % For iteratative and one shot KFE
bayyydelta(:,:,1)     = baydelta;
bayyydelta(:,:,2)     = baydelta;
grid_diag             = spdiags(bayyydelta(:),0,MMMM,MMMM);                % For one shot KFE

% Create mesh: b in 1st dimension, a in 2nd dimension (generally excluded in TS and GT), income type in 2nd/3rd dimension
gridd.bb              = gridd.b*ones(1,NzNed); % repeat b column vector J times in the 2nd dimension
num.Delta             = 1000;                                              % Update parameter in HJB
num.tol_hjb           = 1e-7;                                              % Convergence crietria HJB
etol                  = num.tol_hjb*10;                                    % Relax a little criteria for LCP error (standard practice)
tol_r                 = 1e-7;                                              % Convergence criteria for mkt clearing
maxit                 = 100;                                               % Max number of loops in HJB (outer and inner loop)
num.maxit_hjb         = maxit;
rloopmax              = 40000;                                             % Max number of rloops
lbnd                  = zeros(MMu,1);                                      % LCP bounds
ubnd                  = Inf*ones(MMu,1);                                   % LCP bounds
LCP_errorold          = 1;

%%%%% Build matrix Aswitch summarizing evolution of exogenous states
lparamsC;                                                                  % Pre-allocating memory
born_nodebt;                                                               % Matrix directly flows of newborns (cant be born with debt)
%%%%% Jump matrices
Aswitch1   = [-speye(I)*param.la1,speye(I)*param.la1;speye(I)*param.la2,-speye(I)*param.la2];
Aswitch2   = [-speye(I)*(param.la5+param.la6),speye(I)*param.la5,speye(I)*param.la6;...
             sparse(I,I),-speye(I)*param.la3,speye(I)*param.la3;sparse(I,I),speye(I)*param.la4,-speye(I)*param.la4];
% Initial guess of value functions
lext              = [mean(param.lx3)*ones(I,1);lx1*ones(I,1);lx2*ones(I,1)];
v0pack1           = TSGT_VNONedu(param,guess,num,gridd,Aswitch1,bnodebt);
vne               = [zeros(Ned*I,1),0*repmat(bnodebt*v0pack1.V(:,1),Ned,1)];
v0pack2           = TSGT_Vedu(param,guess,num,gridd,Aswitch2,vne,lext);
v                 = v0pack1.V;
Vstar             = zeros(I,Nz);
cons              = zeros(I,NzNed,abtypes);                                % Pre-allocate memory for V, consumption
Vall              = zeros(I,NzNed,abtypes);
Vallold           = Vall;
i_buy             = min(find(gridd.b+abs(gridd.b(1))>=P0));                % P0 costs i_buy grid points
i_buymesh2        = gridd.b+abs(gridd.b(1))>=P0;
x                 = gridd.b-P0;                                            % Wealth after P0
vold              = ones(I,Nz);
Pold              = P;
buildOmega;
adjMMM2           = zeros(MMM,abtypes);
eEmpeps           = 1e-9;                                                  % For numerical stability - theta_{e,C}>0
param.etol        = etol;
nobdry            = ones(I,Nz);
nobdry(1,:)       = 0;
nobdry            = reshape(nobdry,MMu,1);
nobdryEdu         = ones(I,3);
nobdryEdu(1,:)    = 0;
nobdryEdu         = reshape(nobdryEdu,MMe,1);
ANC               = cell(abtypes,1);                                       % Pre-allocate memory
AC                = cell(abtypes,1);
MMsss             = cell(abtypes,1);
Am                = cell(abtypes,1);   
pointslist        = (1:MMM)';
pointsmat         = reshape(pointslist,I,NzNed);
pointslist0       = (1:I)';
pointsmat0        = reshape(pointslist0, I, 1);
MMsss             = cell(abtypes,1);
Am                = MMsss;
born_rentry       = [sparse(MMM,M-I),repmat(bnodebt,MMM/I,1),sparse(MMM,MMe+M)];
prob_stay         = zeros(MMM,abtypes);
prob_stay(:,1)    = [parent_abtran(2,2)*ones(MMe,1);parent_abtran(4,2)*ones(MMu,1)];
prob_stay(:,2)    = [parent_abtran(1,1)*ones(MMe,1);parent_abtran(3,1)*ones(MMu,1)];
relax             = 0.98;                                                  % Relaxation method parameter
countsmth         = 0;                                                     % Rarely used - See comments in SolveSSEquilibrium
plotsteps         = 100;                                                   % For plotting purposes
comb              = 25000000;                                              % Rarely used - See comments in SolveSSEquilibrium
comb2             = 75;
% Solve for the equilibrium
SolveSSEquilibriumTSGT;

end