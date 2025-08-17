function Res = SSNICLandICLs0(subsidy,GT,Edcostscale,aminfedpctg,SLindic)
% Solves for steady states (SS) of NICL, ICL1 and ICL2
% Gustavo Mellior g.mellior@liverpool.ac.uk

%%%%%%%%% Inputs %%%%%%%%%%
% inpvec(1,1) = Subsidy rate
% inpvec(2,1) = 1 -> GT, 0 -> TS (for SS TSGT for NICL to TSGT without student debt cancellation)
% inpvec(3,1) = For scaling of education costs relative to GDP per capita
% inpvec(4,1) = Baseline student debt limit -> 1
% inpvec(5,1) = Student loan system type: 0 -> ICL2 - 1 -> NICL - 2 -> ICL1

%%
warning('off','all')
Res.codename = mfilename;                                                  % Record file name
fniniguess   = 'startval';                                                 % Baseline parameters mat file
NICL         = SLindic;                                                    % Indicator NICL/ICL1/ICL2
GDPscale     = (328.795/206.538)*54541;                                    % See Appendix
EY           = 0.79;            
P_GDPscale   = (((0.75*21370+0.25*37430)*0.6 + 48510*0.4)/GDPscale)*Edcostscale;
constantkfe  = 0.25;                                                       % Inverse years until graduation
ipc_anddnm   = P_GDPscale*(1/constantkfe);                                 % Edu cost per capita scale
gridd.bmin   = -(20571/GDPscale)*EY;                                       % Overdraft limit
if (aminfedpctg==0)&&(subsidy>0.47)
    %% Load terminal steady state   
    thisScriptPath = fileparts(mfilename('fullpath'));
    NICLandICLsPath = fileparts(thisScriptPath);                           % Go up one level
    mainFolder = fileparts(NICLandICLsPath);                               % Go up again to reach MainFolder
    TSGT_SSres  = fullfile(mainFolder,'TSGT','SSresults');                 % Build the path to TSGT/SSresults
    folderPath  = strcat(TSGT_SSres,'\');
    filepath    = strcat('_bmin_', num2str(abs(round(gridd.bmin, 4)*10000)),'_I_',num2str(205),'_');
    if GT==1
        targetSuffix = strcat(filepath,'GT',num2str(round(100*subsidy)),'GE_','_Edcostscale',num2str(100*round(Edcostscale,2)),'.mat');
    else
        targetSuffix = strcat(filepath,'TS',num2str(round(100*subsidy)),'GE_','_Edcostscale',num2str(100*round(Edcostscale,2)),'.mat');
    end
    SSFile  = load_SSmatching_file(folderPath,targetSuffix);
    curr    = pwd;
    TSGT    = load(fullfile(folderPath,SSFile));
    guess.r = TSGT.guess.r; 
    xi      = TSGT.xi; 
    z_ave   = TSGT.z_ave; 
    KD      = TSGT.KD(end); 
    itax    = TSGT.itax(1); 
    rhoBENCH= TSGT.rho;
    mcBENCH = TSGT.param.mcost(1);
    altr    = TSGT.param.altruism;
    PBENCH  = TSGT.P;
    if GT==1
        itax1 = TSGT.itax(1);
        itax2 = TSGT.itax(2);
    end
    unEmp          = TSGT.unEmp;
    eEmp           = TSGT.eEmp;
    stn            = TSGT.stn;
    rho            = rhoBENCH;
    param.altruism = TSGT.param.altruism;
    param.mcost    = TSGT.mcBENCH*ones(1,2);
else
    load(fniniguess)                                                       % Load baseline parameters
    guess.r = rBENCH;                                                      % Guess interest rate at baseline value
    xi      = xiBENCH;                                                     % CES share
    z_ave   = z_aveBENCH;                                                  % Effective labour supply
    KD      = KBENCH;                                                      % Initial guess capital stock
    itax    = itaxBENCH;                                                   % Initial guess of income tax rate
    if GT==1                                                               % Graduate taxes have two tax rates
        itax1 = itax;                      
        itax2 = itax;
    end
    unEmp          = popBENCH(2);                                          % Initial guess of share in employment with no college education
    eEmp           = popBENCH(5);                                          % Initial guess of share in employment with    college education
    stn            = popBENCH(3);                                          % Initial guess of share in university (students)
    param.altruism = altr;                                                 % Altruism strength
    param.mcost    = mcBENCH*ones(1,2);                                    % Mental cost
    if NICL==1                                                             % Subjective discount rate
        rho        = rhoBENCH*1.02;
    else
        rho        = rhoBENCH;
    end
    if NICL==0
        % As ICL2 takes a very long time to solve, we use a good initial guess
        unEmp          = 0.458940004769962;
        eEmp           = 0.409119798868063;
        stn            = 0.078549550510669;
        itax           = 0.10379;
        KD             = 2.432;
        z_ave          = 0.441;
        guess.r        = 0.03242;
    end
end

param.rho      = rho;
negcons        = 0;                                                        % Track if consumption is negative
%% PARAMETERS AND INITIAL OBJECTS:
if NICL==1
    if (aminfedpctg==0)&&(subsidy>0.47)
        cwp_enforced  = 0;
    else
        cwp_enforced  = 1;                                                 % Enforce CWP=1.7 (for baseline NICL)
    end
else
    cwp_enforced      = 0;
end
update_P_steps        = 5;                                                 % Number of initial rloop iterations after which P is updated
crra                  = 2;                                                 % Risk aversion in CRRA utility
param.s               = crra;                                              % Store in structure used by functions
param.alpha           = 0.36;                                              % Elasticity of Y wrt K
constant              = 0;                                                 % Constant in CRRA utility function
param.constant        = constant;                                          % Store in structure used by functions
constantkfe2          = 1000000;                                           % For KFE jumps (old method - unused)
r                     = guess.r;
onerun                = tic;                                               % Track runtime
param.kappa           = 1/49;                                              % Death rate
subsidy               = 0.47;                                              % Education cost % paid by government
if SLindic==0
    % As ICL2 takes a very long time to solve, we use a good initial guess
    relax       = 0.9999;
    Sloansuplmt = 0.01729;
    P           = 1.304267269039366;
else
    relax       = 0.99;
    P           = PBENCH;
end
Pold                  = P;
SLindic               = NICL;
elastCES              = 2.5;                                               % CES labour
nu                    = (elastCES-1)/elastCES;  
% spread                = 0;
unsubsidized          = 1;                                                 % Unsubsidized loans = 1
if NICL==2
		ataxind       = 1;
		param.amort   = 1/15;                                              % Student loan amortisation rate
		NICL          = 0;
        Picksystem    = 2;
else
		ataxind       = 0;
		param.amort   = 1/15*NICL;                                         % Student loan amortisation rate
        Picksystem    = NICL;
end
lnp                   = (1/30)*(1-NICL);                                   % After 30 years no repayment - Risk premium in student loans (ICL2 only - if NICL then code sets it to zero further below)
studentrepay          = 0;                                                 % 1- Students pay interest and amort , 0- they dont
P0                    = (1-subsidy)*P;                                     % Edu price individuals face
Aprod                 = 1;                                                 % Productivity constant
KYratio               = 3.04;
d                     = param.alpha/KYratio-(0.0505-param.kappa);          % Depreciation
guess.r               = param.alpha/KYratio-d;
r                     = guess.r;
ra                    = guess.r+param.kappa;                               % Interest on student loans - will be updated later, to include cancellation risk
param.ra              = ra;                                                % Store in structure used by functions
% zlow                  = 1;                                                 % Labour efficiency of uneducated
% zhigh                 = zlow*1.7;                                          % Labour efficiency of educated
zstn                  = 0.675*0.75;                                        % Labour efficiency of students
param.ra              = ra; 
unif                  = 0;                    % Zero = nonunif grids on b
param.AssumeConcavity = 0;                    % Do not change - helps for non concavity of V
abtypes               = 2;                    % Ability types
% Transition rates
lx1                   = 0.048;                % Edu skills depreciation rate for unemployed
lx2                   = lx1/2;                % Edu skills depreciation rate for employed
param.lx1             = lx1;                  %
param.lx2             = lx2;                  %
param.lx3_vec         = [1-((0.325+0.519)/2)^(1/4), 1-((0.207+0.256)/2)^(1/4)]; % Strayer and Light (2000) page 315,old and worse than NCES
parent_abtran         = [0.74136275, 0.25863725; 0.34289128, 0.65710872; 0.62180579, 0.37819421; 0.26403233, 0.7359677]; % Intergen trans of ability
param.la1             = 0.157*12;             % Rate from unemp to emp   with no edu 
param.la2             = 0.011*12;             % Rate from emp   to unemp with no edu 
param.la3             = 0.134*12;             % Rate from unemp to emp   with edu 
param.la4             = 0.006*12;             % Rate from emp   to unemp with edu
param.la5             = (1/3)*constantkfe;    % Graduating without a job. Split inverse number of year until graduation into two
param.la6             = (2/3)*constantkfe;    % Graduating with    a job. Split inverse number of year until graduation into two
reprate               = 0.382;                % Unemployment benefits replacement rate
%%%%%%%% Grids %%%%%%%%
% Student loans 
gridd.amin            = -0.7;
gridd.aminfed         = -0.523174541337777*aminfedpctg;
gridd.J               = 10;
J                     = gridd.J;
gridd.amax            = 0; 
gridd.a               = linspace(gridd.amin,gridd.amax,gridd.J); % row vector
gridd.da              = gridd.a(2)-gridd.a(1);  
da                    = gridd.da;
gridd.a               = [fliplr(gridd.amin-gridd.da:-gridd.da:-gridd.da*4+gridd.amin),gridd.a];
gridd.amin            = gridd.a(1);
J                     = numel(gridd.a);
gridd.J               = J;

% Cash
gridd.I               = 204;                 % Number of grid points
g                     = linspace(0.024,0.039,gridd.I);  
new                   = (((1+g).^(1:gridd.I)'-1)./((1+g).^(gridd.I)'-1));
gridd.bmax            = 50;
gridd.b               = [gridd.bmin; gridd.bmin + (gridd.bmax-gridd.bmin).*new];
gridd.I               = numel(gridd.b);      % Number of grid points
I                     = gridd.I;
gridd.bmin            = min(gridd.b(1),0);  
gridd.bzeroind        = max(find(gridd.b<0));% Index for b close to 0

% Income and income taxes
gridd.Nz              = 2;                   % Number of employment states without HE
gridd.Ned             = 3;                   % Number of employment states with HE + student type
Nz                    = gridd.Nz;
NzNed                 = Nz+gridd.Ned;        % Total number of types
gridd.M               = I*J;                 % total number of nodes per type
M                     = gridd.M;
MMu                   = M*Nz;                % Total grid points without HE
gridd.MMu             = MMu;
gridd.MMe             = M*gridd.Ned;
MMe                   = gridd.MMe;           % Total gridd points with HE
MMM                   = MMu+MMe;             % Total grid points
MMMM                  = MMM*abtypes;
% Probability child is born in ability x if parent has edu or not - sss=1
% high ability, sss=2 low ability.
parent_ab_stay        = param.kappa*[[parent_abtran(2,2)*ones(MMe,1);parent_abtran(4,2)*ones(MMu,1)],...
                                     [parent_abtran(1,1)*ones(MMe,1);parent_abtran(3,1)*ones(MMu,1)]]; 
parent_abs_small      = (1-[[parent_abtran(2,2);parent_abtran(4,2)],...
                        [parent_abtran(1,1);parent_abtran(3,1)]]);
parent_ab_stay2       = param.kappa*(1-[[parent_abtran(2,2)*ones(MMe,1);parent_abtran(4,2)*ones(MMu,1)],...
                                        [parent_abtran(1,1)*ones(MMe,1);parent_abtran(3,1)*ones(MMu,1)]]); 
% Ini guess for K demand, wages and taxes
KD                    = (param.alpha*Aprod/(guess.r + d))^(1/(1-param.alpha))*z_ave;  
gridd.wageL           = ((1-param.alpha)*Aprod)*(KD.^param.alpha)*(z_ave^(1-param.alpha-nu))*(xi*(unEmp+zstn*stn)^(nu-1));
gridd.wageH           = ((1-param.alpha)*Aprod)*(KD.^param.alpha)*(z_ave^(1-param.alpha-nu))*((1-xi)*(eEmp)^(nu-1));
gridd.income          = [gridd.wageL*reprate gridd.wageL];          % Unemp and emp no college
gridd.income          = [gridd.income gridd.wageL*zstn gridd.wageH*reprate gridd.wageH]; 
gridd.iiincome        = zeros(gridd.I,gridd.J,NzNed);
param.itax2           = gridd.iiincome;
if GT==1
   param.itax         = [0 itax2 0 0 itax2];     % Guess income tax
else
    param.itax        = [0 itax 0 0 itax];     % Guess income tax
end
for i=1:NzNed
    gridd.iiincome(:,:,i) = gridd.income(1,i);
    param.itax2(:,:,i)    = kron(param.itax(1,i),ones(gridd.I,J));
end

% Delta b (adjustments for nonuniform grids)
dbgrid                = diff(gridd.b);
gridd.dbgridf         = [dbgrid; dbgrid(I-1)];
gridd.dbgridb         = [dbgrid(1); dbgrid];
gridd.dbyF            = repmat(gridd.dbgridf(1:I-1),1,gridd.Ned); 
gridd.dbyB            = repmat(gridd.dbgridb(2:I),1,gridd.Ned);
gridd.dbaF            = repmat(gridd.dbgridf(1:I-1),1,J); 
gridd.dbaB            = repmat(gridd.dbgridb(2:I),1,J);
gridd.dbbayF          = zeros(gridd.I-1,gridd.J,gridd.Ned); 
gridd.dbbayB          = zeros(gridd.I-1,gridd.J,gridd.Ned); 
gridd.dbbayF(:,:,1)   = gridd.dbaF; 
gridd.dbbayF(:,:,2)   = gridd.dbaF;
gridd.dbbayF(:,:,3)   = gridd.dbaF;
gridd.dbbayB(:,:,1)   = gridd.dbaB; 
gridd.dbbayB(:,:,2)   = gridd.dbaB;
gridd.dbbayB(:,:,3)   = gridd.dbaB;
% Trapezoidal rule: for KFE and moments
gridd.bdelta          = zeros(I,1);
gridd.bdelta(1)       = 0.5*dbgrid(1);
gridd.bdelta(2:I-1)   = 0.5*dbgrid(1:I-2) + 0.5*dbgrid(2:I-1);
gridd.bdelta(I)       = 0.5*dbgrid(I-1);
gridd.badelta         = repmat(gridd.bdelta,1,J);
gridd.bydelta         = repmat(gridd.bdelta,1,gridd.Nz);
gridd.baydelta        = zeros(gridd.I,gridd.J,gridd.Nz); 
baydelta              = repmat(gridd.badelta,1,NzNed); % For kdenonuniform3
gridd.baydelta(:,:,1) = gridd.badelta; 
gridd.baydelta(:,:,2) = gridd.badelta;
gridd.bayydelta       = zeros(gridd.I,gridd.J,NzNed); 
gridd.bayydelta(:,:,1:2)  = gridd.baydelta;
gridd.bayydelta(:,:,3:4)  = gridd.baydelta;
gridd.bayydelta(:,:,5)    = gridd.badelta(:,:,1); 
gridd.bayyydelta          = zeros(I,J,NzNed,abtypes);
gridd.bayyydelta(:,:,:,1) = gridd.bayydelta;
gridd.bayyydelta(:,:,:,2) = gridd.bayydelta;
if unif == 1
    gridd.bdelta      = [gridd.bbb(:,1);gridd.bbb(end,1)]; 
    param.aydelta     = repmat(gridd.bdelta,1,gridd.Ned); 
    param.dagridf     = gridd.bdelta;
    param.dagridb     = gridd.bdelta;
    gridd.dbbF        = gridd.bydelta(2:I,:);
    gridd.dbbB        = gridd.bydelta(1:I-1,:);
end
db                    = gridd.bdelta;
% Create mesh: b in 1st dimension, a in 2nd dimension, income type in 3rd dimension
gridd.bb              = gridd.b*ones(1,gridd.J); % repeat b column vector J times in the 2nd dimension
gridd.aa              = ones(gridd.I,1)*gridd.a; % repeat a row vector I times in the 1st dimension
gridd.bbb             = zeros(gridd.I,gridd.J,gridd.Ned); 
gridd.aaa             = zeros(gridd.I,gridd.J,gridd.Ned); 
for hhh=1:3
    gridd.bbb(:,:,hhh) = gridd.bb; 
    gridd.aaa(:,:,hhh) = gridd.aa;
end
la                    = [param.la1,param.la2];% Poisson intensities of non edu groups

% Iteration set-up:
num.Delta             = 1000;                                              % HJB update
num.tol_hjb           = 1e-8;                                              % Convergence crietria HJB
etol                  = num.tol_hjb*100;                                   % Relax a little criteria for LCP error (standard practice)
tol_r                 = 1e-6;                                              % Convergence criteria for mkt clearing
maxit                 = 200;                                               % Max number of loops in HJB (outer and inner loop)
num.maxit_hjb         = maxit;
num.tol_kfe           = 1e-7;                                              % Convergence criteria for KFE
num.mxItkfe           = 100;                                               % Max # of iterations in KFE
num.displaykfe        = 0;                                                 % Set to 1 to see live convergence of KFE
rloopmax              = 5000;
iter                  = 0;
dist                  = zeros(maxit,1);                                    % Store mkt clearing error
d3                    = 1;                                                 % See file printrack and how it updates
d2                    = 1;                                                 % See file printrack and how it updates
lbnd                  = zeros(MMu,1);                                      % LCP bounds
ubnd                  = Inf*ones(MMu,1);                                   % LCP bounds
LCP_errorold          = 1;
% Student loan repayment
param.pctg            = 0.09*(1-NICL);                                     % Tax on earnings in excess of threshold zT (used only in ICL2)
param.zT              = (21000*1.3363)/((66.46/41.253)*37567.76)*EY*(1-NICL); % Earnings threshold (used only in ICL2)
param.atax            = zeros(I,J,NzNed);
bornnodebt_sloans;
% Ini guess for student loan repayment, part of "aDrift"
if ataxind==0 % NICL and ICL2
    param.atax(:,:,1:2)   = -NICL*(param.ra+param.amort)*gridd.aaa(:,:,1:2) + ...
                            (1-NICL)*(((max(gridd.iiincome(:,:,1:2)+guess.r*max(gridd.bbb(:,:,1:2),0)-param.zT,0)))*param.pctg).*(gridd.aaa(:,:,1:2)<0);
    param.atax(:,:,3)     = -studentrepay*NICL*(param.ra+param.amort)*gridd.aaa(:,:,1);        % Students do not pay even if earnings rich
    param.atax(:,:,4)     = -NICL*(param.ra+param.amort)*gridd.aaa(:,:,1) + ...
                            (1-NICL)*(((max(gridd.iiincome(:,:,4)+guess.r*max(gridd.bbb(:,:,1),0)-param.zT,0)))*param.pctg).*(gridd.aaa(:,:,1)<0);
    param.atax(:,:,5)     = -NICL*(param.ra+param.amort)*gridd.aaa(:,:,1) + ...
                            (1-NICL)*(((max(gridd.iiincome(:,:,5)+guess.r*max(gridd.bbb(:,:,2),0)-param.zT,0)))*param.pctg).*(gridd.aaa(:,:,2)<0);
else % ICL1
    param.atax(:,:,1:4)   = zeros(I,J,4);
    param.atax(:,:,5)     = -(param.ra+param.amort)*gridd.aaa(:,:,1);                    	
end                                                               

%% Create matrix capturing flows in student loan space due to graduate tax and employment transitions
tempnp     = [repmat(spdiags(lnp*ones(I,1),0,I,I),J-1,1);sparse(I,I)];       % Auxiliary/temporary variable
lnpvec     = [lnp*ones(I*(J-1),1);zeros(I,1);lnp*ones(I*(J-1),1);zeros(I,1)];% Auxiliary/temporary variable
UKswitchUn = (spdiags(-lnpvec,0,MMu,MMu) +...
             [sparse(M,I*(J-1)),tempnp,sparse(M,M);sparse(M,M),sparse(M,I*(J-1)),tempnp]); % Bswitch for uneducated
lnpvec     = [lnpvec;lnp*ones(I*(J-1),1);zeros(I,1)];
UKswitchEd = (spdiags(-lnpvec,0,MMe,MMe) +...
             [sparse(M,I*(J-1)),tempnp,sparse(M,2*M);sparse(M,M),sparse(M,I*(J-1)),tempnp,sparse(M,M);...
             sparse(M,2*M),sparse(M,I*(J-1)),tempnp]);   % Bswitch for educated                    
Bswitchmaker;
%%
%%%%%% Regions in state space where education is affordable
i_buymesh2   = (gridd.bb-P0 + max(gridd.aa-gridd.aminfed,0))>gridd.b(1);
rloop        = 1;
regionsjump;
bAdj1        = reshape(bAdj(:,:,1),M,1); % What b' does non edu want to keep after paying for edu
bAdj2        = reshape(bAdj(:,:,2),M,1); 
aAdj1        = reshape(aAdj(:,:,1),M,1); % What a' does non edu want to keep after paying for edu
aAdj2        = reshape(aAdj(:,:,2),M,1); 
bAdjvec1     = reshape(bAdj1,M,1);
bAdjvec2     = bAdjvec1;
aAdjvec1     = reshape(aAdj1,M,1);
aAdjvec2     = aAdjvec1;

% Initial guess for value function of non edu
Vall         = zeros(I,J,NzNed,abtypes);
if (aminfedpctg==0)&&(subsidy>0.47)
    for mmm=1:abtypes
        for jjj=1:NzNed
            Vall(:,:,jjj,mmm) = repmat(TSGT.Vall(:,jjj,mmm),1,J,1,1);
        end
    end
    v            = Vall(:,:,1:Nz,1); 
    eEmp_old    = TSGT.eEmp;
    eU_old      = TSGT.eU;
    stn_old     = TSGT.stn;
    unU_old     = TSGT.unU;
    unEmp_old   = TSGT.unEmp;
    stnHab_old  = TSGT.stnHab;
    stnLab_old  = TSGT.stnLab;
else
    V0           = zeros(I,J,Nz);
    V0(:,:,1:2)  = (1-crra)^(-1)*((gridd.iiincome(:,:,1:2).*(1-param.itax2(:,:,1:2))...
                   + (guess.r+param.kappa*(gridd.bbb(:,:,1:Nz)<0)).*gridd.bbb(:,:,1:Nz) - param.atax(:,:,1:2)).^(1-crra)-constant)/(rho+param.kappa);  
    v            = V0;
    guessnonedu;
    v            = V; 
    cons         = zeros(I,J,NzNed,abtypes);
    Vall(:,:,1:Nz,abtypes) = v;
end

vold              = v;
Vbf               = zeros(I,J,Nz);
Vbb               = zeros(I,J,Nz);
vne               = zeros(M*gridd.Ned,1);
Vallold           = Vall;
indx_sss          = zeros(J,Nz,abtypes);
indx_sss1_r       = zeros(J,Nz,rloopmax);
indx_sss2_r       = indx_sss1_r;
grid_diag         = spdiags(gridd.bayyydelta(:)*gridd.da,0,MMMM,MMMM);
adjMMM2           = zeros(MMM,abtypes);
eEmpeps           = 0.01;
adj_r             = zeros(MMu,rloopmax,abtypes);
comb              = rloopmax*0.9;
comb2             = 7;
indxpts           = repmat((1:I)',1,J);
row_r             = zeros(J,Nz,abtypes,rloopmax);
adjmat0           = zeros(I,J,Nz);
gridd.btilde      = max(gridd.bb+gridd.aa,gridd.b(1));
gridd.atilde      = gridd.aa+gridd.bb-gridd.btilde;
bAdjAux           = zeros(I,J,gridd.Nz);
aAdjAux           = zeros(I,J,gridd.Nz);
VstarAux_         = zeros(I,J,gridd.Nz);
param.etol        = etol;
nobdry            = ones(I,J,2);
nobdry(1,:,:)     = 0;
nobdry(:,J,:)     = 0;
nobdry            = reshape(nobdry,MMu,1);
nobdryEdu         = ones(I,J,3);
nobdryEdu(1,:,:)  = 0;
nobdryEdu(:,J,:)  = 0;
nobdryEdu         = reshape(nobdryEdu,MMe,1);
ANC               = cell(abtypes,1);
AC                = cell(abtypes,1);
MMsss             = cell(abtypes,1);
Am                = cell(abtypes,1);   
pointslist        = (1:MMM)';
pointsmat         = reshape(pointslist, I, J, NzNed);
MMsss             = cell(abtypes,1);
Am                = MMsss;
born_rentry       = [sparse(MMM,M-I),repmat(bnodebt,MMM/I,1),sparse(MMM,MMe+M)];
prob_stay         = zeros(MMM,abtypes);
prob_stay(:,1)    = [parent_abtran(2,2)*ones(MMe,1);parent_abtran(4,2)*ones(MMu,1)];
prob_stay(:,2)    = [parent_abtran(1,1)*ones(MMe,1);parent_abtran(3,1)*ones(MMu,1)];
wL                = gridd.wageL;
wH                = gridd.wageH;
dfdt              = 0.0001;
mommatch          = NICL*cwp_enforced;
vdistool          = 1;
pltsteps          = 10;

deldec1_twoassets_TSGT_2025_ICL2testAug9b
end