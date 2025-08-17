
for rloop=1:rloopmax      % Interest rate loop
    r_r(rloop) = guess.r;
    okooloop = 0;
    movens   = zeros(abtypes,1);
    switch NICL
        case 1
            try
                guess    = rmfield(guess,'vprev');
            catch
            end
            % Perturbing the initial guess of the value function at every
            % rloop helps in NICL. Alternatively increase the relax parameter
            Vall     = Vall + 0.001*repmat(gridd.bb,1,1,NzNed,abtypes) + 0.001*repmat(gridd.aa,1,1,NzNed,abtypes);
    otherwise
    end
for ooloop=1:15
    
    for sss=1:abtypes
            param.lx3                 = param.lx3_vec(sss);
            computeV_sl;
            if (okooloop==1)||(ooloop>14)
                cons(:,:,1:Nz,sss)       = c;                                  % Store consumption policy
                cons(:,:,Nz+1:NzNed,sss) = v0pack2.c;
                ERBetterNonedu        = abs(VstarER(:)-V_stacked)<param.etol;  % Early repayment better than staying put
                UniBetterER           = Vstar_GoToUni>VstarER(:);              % University better than early repayment
                UniBetterNoedu        = Vstar_GoToUni>V_stacked;               % University better than staying put
                adj                   = abs(V_stacked - Vstar_GoToUni)<param.etol;      % Going to university - Allow for some tolerance in defining solution for adjustment decision
                adjall                = (abs(V_stacked - VstarAll_stacked)<param.etol); % Any adjustment - Allow for some tolerance in defining solution for adjustment decision
                adjmat                = reshape(adj,I,J,Nz);
                adjmat2               = zeros(I,J,Nz);
                adjmat0(:,:,1)        = (adjmat(:,:,1)*1).*indxpts;
                adjmat0(:,:,2)        = (adjmat(:,:,2)*1).*indxpts;
                adjmat0(adjmat0==0)   = I+1;
                % On the upper end of the "go to university region" a sparse
                % grid may cycle between neighbouring points as we get close to equilibrium. As these regions
                % have little mass, averaging out the grid boundaries reduces computational costs. If rloop>comb, then we avoid cycling 
                % by averaging (rarely used). Alternatively, increase grid resolution. 
                if rloop>comb
                    adj               = round((1/comb2)*( adj + sum( adj_r(:,rloop-comb2+1:rloop-1,sss),2)  ));
                end
                adj                   = logical(adj.*repmat(i_buymesh2(:),2,1));
                adjMMM                = [adj;false(MMe,1)];
                adjRegion1            = reshape(adj(1:M),I,J).*ones(I,J);        % Transforms logicals into floats
                adjRegion2            = reshape(adj(M+1:MMu),I,J).*ones(I,J);
                adjRegion(:,:,1)      = adjRegion1; 
                adjRegion(:,:,2)      = adjRegion2; 
                adjMMM2(:,sss)        = adjMMM;
                adjMMM2_ER(:,sss)     = [ERBetterNonedu.*~adj.*nobdry;v0pack2.RepayEarly.*nobdryEdu];
                adjMMM2_ERb(:,sss)    = [ERBetterNonedu.*adjall.*nobdry;v0pack2.RepayEarly.*nobdryEdu];
                adj_r(:,rloop,sss)    = adj; 
                ANC{sss}              = A;
                AC{sss}               = v0pack2.Amat;                      % Store parts of A matrix (by ability type)
                movens(sss)           = 1;
            end
        end
    Vallnew          = Vall;
    vdistool(ooloop) = max(abs(Vallnew(:)-Vallold(:)));
    Vallold          = Vallnew;
    d4(:,sss,rloop)  = [d2;vdistool(ooloop)];
    clc
    if vdistool(ooloop)<num.tol_hjb
        okooloop = 1;
        if sum(movens)>1
            break
        end
    end

end

getdistribution_NICLandICL;

if max(abs(sum(Afinal,2)))>1e-9
    figure(1)
    plot(sum(Afinal,2))
    grid on
    drawnow
    disp('Check flows in A Matrix')
end

unU          = gridd.da*gridd.bdelta'*sum(sum(g(:,:,1,:),4),2);    % Share unemployed with no education
unEmp        = gridd.da*gridd.bdelta'*sum(sum(g(:,:,2,:),4),2);    % Share employed with no education
stn          = gridd.da*gridd.bdelta'*sum(sum(g(:,:,3,:),4),2);    % Share student
eU           = gridd.da*gridd.bdelta'*sum(sum(g(:,:,4,:),4),2);    % Share unemployed WITH education
eEmp         = gridd.da*gridd.bdelta'*sum(sum(g(:,:,5,:),4),2);    % Share employed WITH education
stnHab       = gridd.da*gridd.bdelta'*sum(g(:,:,3,1),2);           % Share student high ability 
stnLab       = gridd.da*gridd.bdelta'*sum(g(:,:,3,2),2);           % Share student low ability 
t5Hab        = gridd.da*gridd.bdelta'*sum(sum(g(:,:,[4 5],1),3),2);% Share graduates high ability 
t5Lab        = gridd.da*gridd.bdelta'*sum(sum(g(:,:,[4 5],2),3),2);% Share graduates low ability 

switch NICL 
    case 1 % Enforce CWP at baseline
         xi   = ((max(eEmp,eEmpeps)/(unEmp+zstn*stn))^(nu-1))/(cwp+(max(eEmp,eEmpeps)/(unEmp+zstn*stn))^(nu-1)); % CES share
    otherwise
end

z_ave        = ( xi*((unEmp+zstn*stn)^nu) + (1-xi)*((max(eEmp,eEmpeps))^nu) )^(1/nu); % Effective labour supply
studentdebt  = gridd.da*(sum(sum(sum(g.*gridd.bayyydelta,4),3),1))*gridd.a';          % Aggregate student debt
KS(rloop)    = gridd.da*(gridd.b.*gridd.bdelta)'*sum(sum(sum(g,4),3),2) + studentdebt;% Capital supply 
rnew         = (z_ave.^(1-param.alpha)).*(param.alpha* Aprod*KS(rloop).^(param.alpha-1)) - d; % Interest rate
wL_n         = min(((1-param.alpha)*Aprod)*((KS(rloop).^param.alpha).*(z_ave.^(1-param.alpha-nu))).*(xi*(unEmp+zstn*stn).^(nu-1)),1); % Wages
wH_n         = min(max(((1-param.alpha)*Aprod)*((KS(rloop).^param.alpha).*(z_ave.^(1-param.alpha-nu))).*((1-xi)*(max(eEmp,eEmpeps)).^(nu-1)),wL_n),2*wL_n);
KD(rloop)    = (param.alpha*Aprod/(guess.r + d))^(1/(1-param.alpha))*(z_ave); % Capital demand
SS(rloop)    = KS(rloop) - KD(rloop);
switch NICL
    case 0
        r         = relax.*r + (1-relax).*rnew;
        guess.r   = r;
    case 1
        if abs(SS(rloop))<1
            rho   = min(max(rho + 10*dfdt*SS(rloop),0.0001),0.095);
        else
            rho   = min(max(rho + 0.001*dfdt*SS(rloop),0.0001),0.095);
        end
        param.rho = rho;
        guess.r   = (param.alpha/KYratio-d); % Consistent with baseline K/Y ratio
        r         = guess.r;
end
wL           = relax.*wL + (1-relax).*wL_n;
wH           = relax.*wH + (1-relax).*wH_n;
gridd.wageL  = wL;
gridd.wageH  = wH;
AggC         = gridd.da*gridd.bdelta'*(sum(sum(sum(g.*cons,4),3),2));      % Aggregate consumption
KD(rloop)    = (param.alpha*Aprod/(guess.r + d))^(1/(1-param.alpha))*(z_ave); % Capital demand
AggY         = Aprod*(KD(rloop)^param.alpha)*(z_ave)^(1-param.alpha);      % Aggregate output
ra           = guess.r+param.kappa; 
param.ra     = ra; 
% Update income grid
gridd.income = [gridd.wageL*reprate gridd.wageL];                          
gridd.income = [gridd.income gridd.wageL*zstn gridd.wageH*reprate gridd.wageH]; 
for jjj=1:NzNed
    gridd.iiincome(:,:,jjj) = gridd.income(1,jjj)*ones(gridd.I,gridd.J,1);
end
% Update student loan tax
if ataxind==0
    param.atax(:,:,1:2)   = -NICL*(param.ra+param.amort)*gridd.aaa(:,:,1:2) + ...
                            (1-NICL)*(((max(gridd.iiincome(:,:,1:2)+guess.r*max(gridd.bbb(:,:,1:2),0)-param.zT,0)))*param.pctg).*(gridd.aaa(:,:,1:2)<0);
    param.atax(:,:,3)     = -studentrepay*NICL*(param.ra+param.amort)*gridd.aaa(:,:,1); 
    param.atax(:,:,4)     = -NICL*(param.ra+param.amort)*gridd.aaa(:,:,1) + ...
                            (1-NICL)*(((max(gridd.iiincome(:,:,4)+guess.r*max(gridd.bbb(:,:,1),0)-param.zT,0)))*param.pctg).*(gridd.aaa(:,:,1)<0);
    param.atax(:,:,5)     = -NICL*(param.ra+param.amort)*gridd.aaa(:,:,1) + ...
                            (1-NICL)*(((max(gridd.iiincome(:,:,5)+guess.r*max(gridd.bbb(:,:,2),0)-param.zT,0)))*param.pctg).*(gridd.aaa(:,:,2)<0);
else
    param.atax(:,:,1:4)   = 0*repmat(param.amort*gridd.aa,1,1,4);
    param.atax(:,:,5)     = -(param.ra+param.amort)*gridd.aaa(:,:,1);                        
end
% Cost of student loan system
switch Picksystem
    case 0
        Sloansuplmt0 = -ra*gridd.da*(sum(sum(sum(g(:,:,[1,2,4],:).*gridd.bayyydelta(:,:,[1 2 4],:),4),3),1))*gridd.a';
        Sloansuplmt1 = ra*gridd.da*gridd.aa(:,1)'*(sum(g(:,1,3,1).*gridd.bdelta(:,1),2)) + ra*gridd.da*gridd.aa(:,1)'*(sum(g(:,1,5,1).*gridd.bdelta(:,1),2))...
                      +ra*gridd.da*gridd.aa(:,1)'*(sum(g(:,1,3,2).*gridd.bdelta(:,1),2)) + ra*gridd.da*gridd.aa(:,1)'*(sum(g(:,1,5,2).*gridd.bdelta(:,1),2)); % interest subsidy at \bar{\bar{a}}
        Sloansuplmtnew = Sloansuplmt0 - Sloansuplmt1 + lnp*(-studentdebt); 
        Sloansuplmt    = Sloansuplmtnew*(1-relax) + Sloansuplmt*relax;
    case 1
        Sloansuplmt  = unsubsidized*(-ra*gridd.da*gridd.aa(:,1)'*( sum(g(:,1,3,:),4).*gridd.bayydelta(:,1,3).*(aDrift(:,1,3)<0) ))...
                      + (1-unsubsidized)*(-ra*gridd.da*( sum(sum(sum(repmat(gridd.aa,1,1,1,abtypes).*g(:,:,3,:).*gridd.bayyydelta(:,:,3,:),4),2),1) ));
    case 2
        Sloansuplmt  = - sum((ra+param.amort)*gridd.da*( sum(sum(sum(repmat(gridd.aa,1,1,NzNed-1,abtypes).*g(:,:,[1:4],:).*gridd.bayyydelta(:,:,3,:),4),2),1) ))...
                       + lnp*(-studentdebt);
end

Bswitchmaker;

if min(gridd.iiincome(1,1,:) + guess.r*gridd.bbb(1,1,1) - param.atax(1,1,:)) < 0
    % Track if debt limit is too lose
    negcons = negcons + 1;   
end
if negcons>5
    % If the lowest consumption is negative too many times, exit, change debt limits
    break
end

% Update education costs and regions where education is affordable
if rloop>update_P_steps
    switch SLindic
        case 0
        otherwise
            P    = 0.5*(ipc_anddnm*(AggY)) + 0.5*Pold;  
    end
    Pold         = P;
    P0           = (1-subsidy)*P;                                          % Education price individuals face
    i_buymesh2   = (gridd.bb -P0 + max(gridd.aa-gridd.aminfed,0))>gridd.b(1);
    regionsjump;
    bAdj1       = reshape(bAdj(:,:,1),M,1);                                % What b' does non edu want to keep after paying for edu
    bAdj2       = reshape(bAdj(:,:,2),M,1); 
    aAdj1       = reshape(aAdj(:,:,1),M,1);                                % What a' does non edu want to keep after paying for edu
    aAdj2       = reshape(aAdj(:,:,2),M,1); 
    bAdjvec1    = reshape(bAdj1,M,1);       
    bAdjvec2    = bAdjvec1;                 
    aAdjvec1    = reshape(aAdj1,M,1);       
    aAdjvec2    = aAdjvec1;                 
end

% Education costs and income tax
unbenef      = reprate*(unU*gridd.wageL + eU*gridd.wageH);                 % Unemployment benefits
flowsINEDUn1 = param.lx1*eU + param.lx2*eEmp + param.lx3_vec(1)*stnHab + param.lx3_vec(2)*stnLab;
flowsINEDUn2 = (eU + eEmp + stn)*param.kappa;
edinflow_old = flowsINEDUn1+flowsINEDUn2;
edcost       = (flowsINEDUn1+flowsINEDUn2)*P;                              % Education costs
subs         = subsidy*edcost;                                             % Public education costs
MKT          = AggY - AggC - d*KS(rloop) - edcost;                         % Resource const error
switch GT
    case 0
        itaxnew    = (unbenef + subs + Sloansuplmt)/(gridd.wageL*unEmp + gridd.wageH*eEmp);
        itax       = relax*itax + (1-relax)*itaxnew;
        itax       = min(0.5,itax);
        param.itax = [0 itax 0 0 itax];
    case 1
        itaxnew1   = unbenef/(gridd.wageL*unEmp + gridd.wageH*eEmp);
        itaxnew2   = itaxnew1 + ((subs + Sloansuplmt)/(gridd.wageH*eEmp))*(eEmp>=eEmpeps);
        itax1      = relax*itax1 + (1-relax)*itaxnew1;
        itax2      = relax*itax2 + (1-relax)*itaxnew2;
        itax1      = min(0.5,itax1);
        itax2      = min(0.5,itax2);
        param.itax = [0 itax1 0 0 itax2];
end

for jjj=1:NzNed
    param.itax2(:,:,jjj) = kron(param.itax(1,jjj),ones(gridd.I,J));
end

% Track evolution of results
eEmp_r(rloop)    = eEmp+eU;
P_r(rloop)       = P;
w_r(rloop,:)     = [gridd.wageL,gridd.wageH];
Y_r(rloop,1)     = AggY;
share_r(rloop,:) = [unU,unEmp,eU,eEmp];
MKT_r(rloop,1)   = MKT;
C_r(rloop,:)     = AggC;
itax_r(rloop,:)  = param.itax(end);
r_r(rloop,:)     = guess.r;
z_ave_r(:,rloop) = z_ave;

if (abs(SS(rloop))<tol_r)&&(rloop>15)&&(abs(MKT_r(end))<tol_r) 
    break
end

end
Res.Welfare   = gridd.da*sum(sum(sum(Vall.*g,4),3),2)'*gridd.bdelta;
inequalitystats;
Res.K         = KD(end);
Res.negcons   = negcons;
Res.r         = guess.r;
Res.SSconv    = SS(end);
Res.mktclr    = MKT_r(end);
Res.convNed   = d1;
Res.convint   = d2;
Res.Y         = AggY;
Res.tax       = itax;
Res.rloop     = rloop;
Res.oloop     = oloop;
Res.iloopNed  = n;
Res.Meanc     = AggC;
Res.Varc      = Varc;
Res.Giniw     = Giniw;
Res.Giniwa    = Giniw2;
Res.Ginic     = Ginic;
Res.CVc       = CVc;
Res.CVw       = CVb;
Res.CVa       = CVa;
Res.timec     = toc(onerun);
Res.debtors   = Bdebtors;
Res.adebtors  = Adebtors;
Res.abdebtors = ABdebtors;
Res.pop       = [unU,unEmp,stn,eU,eEmp];
Res.abtypes   = sum(gridd.da*repmat(gridd.bdelta,1,2).*squeeze(sum(sum(g,3),2)),1);
Res.stnabpop  = [stnHab,stnLab];
Res.gradabpop = [t5Hab,t5Lab];
Res.wage      = [gridd.wageL,gridd.wageH];
Res.mnlinc    = Res.pop*gridd.income';                                     
Res.netearn   = Res.pop*(gridd.income.*(1-itax))';
Res.UItax     = unbenef/(gridd.wageL*unEmp + gridd.wageH*eEmp);
Res.z_ave     = z_ave;
Res.VV        = Vall;
Res.Cbig      = AfT';
Res.dens      = g;
Res.date      = clock;
nus           = ceil(nu*1000)/100;
Res.benchmark_NICLparams = fniniguess;

% Create 'SSresults' folder if it doesn't exist
foldername = fullfile(cd,'SSresults');
if ~exist(foldername,'dir')
    mkdir(foldername);
end

currentDate = datestr(now,'mmmddyyyy'); 
prefname    = fullfile(foldername,currentDate);
filepath    = strcat(prefname,'_bmin_',num2str(abs(round(gridd.b(1),4)*10000)),'_I_',num2str(gridd.I),'_amin_',num2str(abs(round(gridd.a(1),3)*1000)),'_J_',num2str(J));

if aminfedpctg==0
    if GT==0
        matfilenm = strcat(filepath,'TS',num2str(round(100*subsidy)),'GE_','_Edcostscale',num2str(100*round(Edcostscale,2)));
    else
        matfilenm = strcat(filepath,'GT',num2str(round(100*subsidy)),'GE_','_Edcostscale',num2str(100*round(Edcostscale,2)));
    end
else
    if SLindic==0
        matfilenm = strcat(filepath,'_ICL2GE_','Edcostscale',num2str(100*round(Edcostscale,2)),'_aminfedpctg',num2str(100*aminfedpctg),'_amort',num2str(round(1/param.amort)),'_subs',num2str(100*subsidy));
    elseif SLindic==1
        matfilenm = strcat(filepath,'_NICLGE_','Edcostscale',num2str(100*round(Edcostscale,2)),'_aminfedpctg',num2str(100*aminfedpctg),'_amort',num2str(round(1/param.amort)),'_subs',num2str(100*subsidy));
    else
        matfilenm = strcat(filepath,'_ICL1GE_','Edcostscale',num2str(100*round(Edcostscale,2)),'_aminfedpctg',num2str(100*aminfedpctg),'_amort',num2str(round(1/param.amort)),'_subs',num2str(100*subsidy));
    end
end
save(matfilenm)

% For moment matching, save baseline values
% if NICL==1
% % % Save Benchmark values (GE NICL)
%     xiBENCH    = xi;
%     z_aveBENCH = z_ave;
%     rBENCH     = r;
%     popBENCH   = Res.pop;
%     itaxBENCH  = itax;
%     PBENCH     = P;
%     KBENCH     = KD(end);
%     rhoBENCH   = rho;
%     mcBENCH    = param.mcost(1);
%     altr       = param.altruism;
%     dBENCH     = d;
%     save startval xiBENCH z_aveBENCH rBENCH popBENCH itaxBENCH PBENCH KBENCH rhoBENCH mcBENCH altr dBENCH
% end