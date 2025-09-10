function Resret = NICLtoTSGT(inpts,GTindic,gpoints)

filenametran = mfilename;    % Record the matlab file that produces results
%% Load initial steady state
load('../SteadyStates/NICLandICLs/SSresults/Aug172025_bmin_1872_I_205_amin_1011_J_14_NICLGE_Edcostscale100_aminfedpctg100_amort15_subs47');
gg0      = reshape(squeeze(sum(g*gridd.da,2)),I*NzNed*abtypes,1);
r00      = guess.r;                                                        %#ok
K0       = KD(end);                                                        %#ok
W00      = Res.Welfare;                                                    %#ok
xi0      = xi;                                                             %#ok
Pini     = P;                                                              %#ok
amort    = param.amort;                                                    %#ok
itax0    = itax;                                                           %#ok
unEmp0   = unEmp;                                                          %#ok
stn0     = stn;                                                            %#ok
eEmp0    = eEmp;                                                           %#ok
unU0     = unU;                                                            %#ok
Giniw0   = Res.Giniwa;                                                     %#ok
CVw0     = Res.CVw;                                                        %#ok
meaninc0 = Res.mnlinc;                                                     %#ok
debtors0 = Res.abdebtors;                                                  %#ok
wage0    = Res.wage;                                                       %#ok
origI    = I;
origb    = gridd.b;

%% Load terminal steady state   
folderPath = '../SteadyStates/TSGT/SSresults/';
if GTindic==0
    targetSuffix = strcat('_bmin_1872_I_',num2str(gpoints),'_TS',num2str(inpts(1)*100),'GE__Edcostscale100.mat');
else
    targetSuffix = strcat('_bmin_1872_I_',num2str(gpoints),'_GT',num2str(inpts(1)*100),'GE__Edcostscale100.mat');
end
SSFile  = load_SSmatching_file(folderPath,targetSuffix);
load(fullfile(folderPath,SSFile));
addpath('../SteadyStates/TSGT/Functions');
clear('diff');clear('adjmat0');
if I~=205
    gg0orig             = gg0;
    gorig               = reshape(gg0orig,origI,NzNed,abtypes);
    for iii=1:NzNed
        for jjj=1:abtypes
            g0interp(:,iii,jjj) = interp1(origb,gorig(:,iii,jjj),gridd.b);
        end
    end
    gg0                 = reshape(g0interp./sum(g0interp.*bayyydelta,'all'),MMMM,1);
    gg0                 = sparse(gg0);
end
edinflow            = (flowsINEDUn1+flowsINEDUn2);
edinflowT           = edinflow;                                            %#ok
K_st                = KD(end);                                             %#ok
gT                  = g;                                            
v_st                = Vall;                                                %#ok
titaxend            = itax;                                                %#ok
tg                  = 0.04;
Time                = 0.001 + (640-0.001).*(((1+tg).^(1:185)'-1)/((1+tg).^185-1));
dtsvec              = diff(Time);                                          
T                   = Time(end);                                           %#ok
perds               = length(Time);                                        % Time grid points
dtsvec              = [dtsvec(1);dtsvec];                                  %#ok
timeref;                                                                   % Refine time grid
r_t                 = guess.r*ones(perds,1);
guess.r_t           = r_t;
num.maxit_hjb       = 1;
param.looplimit     = 1;
z_ave_t             = ones(perds,1)*z_ave;
z_ave_old           = z_ave_t;                                             %#ok
relax               = 0.995*ones(perds,1);                                 % Path update parameter
K_t                 = K_st*ones(perds,1);                                  %Initial guess of the path of K_t
itax_t              = kron(itax',ones(perds,1));                           
itax_told           = itax_t;                                              %#ok
V_t                 = zeros(I,NzNed,abtypes,perds);                        %#ok
cT                  = cons;                                                %#ok
laction             = cell(perds,1);                                       %#ok
AggKdot_t           = zeros(perds,1);                                      
z_ave_T             = z_ave;                                               %#ok
eEmp_t              = eEmp*ones(perds,1);                                  %#ok
stn_t               = stn*ones(perds,1);                                   %#ok
unU_t               = unU*ones(perds,1);                                   %#ok
unEmp_t             = unEmp*ones(perds,1);                                 %#ok
eU_t                = eU*ones(perds,1);                                    %#ok
stnHab_t            = (gridd.bdelta'*g(:,3,1))*ones(perds,1);              %#ok
stnLab_t            = (gridd.bdelta'*g(:,3,2))*ones(perds,1);              %#ok           
unU_T               = unU;                                                 %#ok
unEmp_T             = unEmp;                                               %#ok
stn_T               = stn;                                                 %#ok
eU_T                = eU;                                                  %#ok
eEmp_T              = eEmp;                                                %#ok
Giniw_t             = NaN(perds,1);                                        %#ok
Ginic_t             = NaN(perds,1);                                        %#ok
debtors_t           = NaN(perds,1);                                        %#ok
AggC_t              = AggKdot_t;                                           %#ok
AggY_t              = AggKdot_t;                                           %#ok
CVc_t               = AggKdot_t;                                           %#ok
CVw_t               = AggKdot_t;                                           %#ok
cons_t              = zeros(I,NzNed,abtypes,perds);                        %#ok
P_t                 = P*ones(perds,1);
P_told              = P_t;                                                 %#ok
meaninc_t           = AggKdot_t;                                           %#ok
maxit               = 1;                                                   %#ok
num.maxit_hjb       = 1;                                                   %#ok
C_t                 = cell(perds,1);                                       %#ok
bayyydeltavec       = bayyydelta(:);                                       %#ok
Aft                 = Afinal;                                              %#ok
nu                  = (elastCES-1)/elastCES;                               %#ok
iniedjump           = zeros(perds,1);                                      %#ok
Vallold             = ones(MMMM,1);                                        %#ok
loopcount           = zeros(3,2,perds);                                    %#ok
indxpts             = (1:I)';                                              %#ok
comb3               = 1;                                                   %#ok
comb4               = 6;                                                   %#ok
adj_tr              = zeros(MMu,comb4,abtypes,perds);                      %#ok - For stability, see steady state files
for sss=1:abtypes
    adj_tr(:,:,sss,:) = repmat(adjMMM2(1:MMu,sss).*repmat(i_buymesh2(:),2,1),1,comb4,1,perds);
end
pathA               = studentdebt*exp(-amort*Time);                        %#ok
probHHnc            = 1-parent_abs_small(1,1);                             %#ok
probHHc             = 1-parent_abs_small(2,1);                             %#ok
probHLnc            = parent_abs_small(1,1);                               %#ok
probHLc             = parent_abs_small(2,1);                               %#ok
probLLnc            = 1-parent_abs_small(1,2);                             %#ok
probLLc             = 1-parent_abs_small(2,2);                             %#ok
probLHnc            = parent_abs_small(1,2);                               %#ok
probLHc             = parent_abs_small(2,2);                               %#ok
wL_t                = gridd.wageL*ones(perds,1);                           %#ok
wH_t                = gridd.wageH*ones(perds,1);                           %#ok
edinflow_t_old      = edinflowT*ones(perds,1);                             %#ok
edinflow_t          = edinflow_t_old;                                      %#ok
clear('gg')
clear('omegagen')

%% COMPUTE THE TRANSITION
SolveNICLtoTSGT;

%% Save results
% Create 'SSresults' folder if it doesn't exist
foldername = fullfile(cd,'Transitionresults');
if ~exist(foldername,'dir')
    mkdir(foldername);
end

currentDate = datestr(now,'mmmddyyyy'); 
dateinname  = fullfile(foldername,currentDate);

if GTindic==1
    transfname = strcat(dateinname,'_NICL_to_GT_bmin_',num2str(abs(round(-gridd.b(1),4)*10000)),'_subs_',num2str(inpts(1)*100));
else
    transfname = strcat(dateinname,'_NICL_to_TS_bmin_',num2str(abs(round(-gridd.b(1),4)*10000)),'_subs_',num2str(inpts(1)*100));
end
save(transfname)
Resret       = checkconv_r(rloop);
end
