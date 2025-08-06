
% Prices loop
for rloop=1:rloopmax
    r_r(rloop)     = guess.r;
    % Perturb the previous solution to avoid numerical errors with LCP -
    % Use new initial guess for each price loop
    clear('guess.vprev')
    v0pack1        = TSGT_VNONedu(param,guess,num,gridd,Aswitch1,bnodebt);
    v              = v0pack1.V;
    Vall(:,1:Nz,1) = v;
    Vall(:,1:Nz,2) = v;
    okooloop       = 0; 
    % HJBVI loop
for ooloop=1:10
    
    for sss=1:abtypes
        param.mcost = mcost(sss);
        param.lx3   = param.lx3_vec(sss);
        lext        = kron([param.lx3;param.lx1;param.lx2],ones(M,1));
        computeVs_TSGT_for_ooloop;
        if (okooloop==1)||(ooloop>9)
            cons(:,1:Nz,sss)       = c;
            cons(:,Nz+1:NzNed,sss) = v0pack2.c;
            adj                    = abs(Vstar_stacked-V_stacked)<etol;    % Region where agents enrol in university
            adjmat                 = reshape(adj,I,Nz);
            % On the upper end of the "go to university region" a sparse
            % grid may cycle between neighbouring points as we get close to equilibrium. As these regions
            % have little mass, averaging out the grid boundaries reduces computational costs. If rloop>comb, then we avoid cycling 
            % by averaging (rarely used). Alternatively, increase grid resolution. 
            switch smallgrid
                case 1
                    if rloop>comb
                        for jjj=1:J
                            adj           = round((1/comb2)*( adj + sum( adj_r(:,rloop-comb2+1:rloop-1,sss),2)  ));
                        end
                    end
            end
            adj                      = logical(adj.*repmat(i_buymesh2(:),2,1)); % Added to avoid choosing edu when not affordable
            adjMMM                   = [adj;false(MMe,1)];
            adjRegion1               = reshape(adj(1:M),I,J).*ones(I,J);        % Multiplication by ones(I,J) transforms logicals into floats
            adjRegion2               = reshape(adj(M+1:MMu),I,J).*ones(I,J);
            adjRegion(:,1)           = adjRegion1; 
            adjRegion(:,2)           = adjRegion2; 
            adjMMM2(:,sss)           = adjMMM;
            adj_r(:,rloop,sss)       = adj; 
            ANC{sss}                 = A;                                  % A matrix no edu
            AC{sss}                  = v0pack2.Amat;                       % A matrix students and edu
            Vstar_stacked_sss(:,sss) = Vstar_stacked;
        end
        
    end
    Vallnew          = Vall;
    vdistool(ooloop) = max(abs(Vallnew(:)-Vallold(:)));
    if (okooloop==1)
        break
    end
    if vdistool(ooloop)<num.tol_hjb
        okooloop = 1;
%         disp('****************** OLOOP CONVERGED ****************')
%         break
    end
    Vallold = Vallnew;
end

getdistribution_SSTSGT;
unU        = gridd.bdelta'*sum(g(:,1,:),3);                                % Share unemployed with no education
unEmp      = gridd.bdelta'*sum(g(:,2,:),3);                                % Share employed with no education
stn        = gridd.bdelta'*sum(g(:,3,:),3);                                % Share student
eU         = gridd.bdelta'*sum(g(:,4,:),3);                                % Share unemployed WITH education
eEmp       = gridd.bdelta'*sum(g(:,5,:),3);                                % Share employed WITH education
stnHab     = gridd.bdelta'*g(:,3,1);                                       % Share student high ability (low dropout)
stnLab     = gridd.bdelta'*g(:,3,2);                                       % Share student low ability 
z_ave        = ( xi*((unEmp+zstn*stn)^nu) + (1-xi)*((max(eEmp,eEmpeps))^nu) )^(1/nu); % Effective labour supply
KS(rloop)    = max((gridd.bdelta.*gridd.b)'*sum(sum(g,2),3),0.5);          % Capital supply 
rnew         = (z_ave.^(1-param.alpha)).*(param.alpha* Aprod*KS(rloop).^(param.alpha-1)) - d; % Interest rates
wL_n         = min(((1-param.alpha)*Aprod)*((KS(rloop).^param.alpha).*(z_ave.^(1-param.alpha-nu))).*(xi*(unEmp+zstn*stn).^(nu-1)),1); % Wages
wH_n         = min(max(((1-param.alpha)*Aprod)*((KS(rloop).^param.alpha).*(z_ave.^(1-param.alpha-nu))).*((1-xi)*(max(eEmp,eEmpeps)).^(nu-1)),wL_n),2*wL_n);
Kdiff(rloop) = abs(KS(rloop)-Kold);
Kold         = relax.*Kold  + (1-relax).*KS(rloop);                        % Slow update of K, r, w
r            = relax.*r     + (1-relax).*rnew;
guess.r      = r;
wL           = relax.*wL + (1-relax).*wL_n;
wH           = relax.*wH + (1-relax).*wH_n;
KD(rloop)    = (param.alpha*Aprod/(guess.r + d))^(1/(1-param.alpha))*(z_ave); % Capital demand
SS(rloop)    = KS(rloop) - KD(rloop);
AggC         = gridd.bdelta'*(sum(sum(g.*cons,3),2));
AggY         = Aprod*(KD(rloop)^param.alpha)*(z_ave)^(1-param.alpha);

% Update wage
gridd.wageL  = wL;
gridd.wageH  = wH; 
gridd.income = [gridd.wageL*reprate,gridd.wageL,gridd.wageL*zstn,gridd.wageH*reprate,gridd.wageH]; 
for jjj=1:NzNed
    gridd.iiincome(:,jjj) = gridd.income(1,jjj)*ones(gridd.I,1);
end

if min(gridd.iiincome(1,1,:) + (guess.r+param.kappa)*gridd.bb(1,1,1)) < 0
    disp('WARNING: borrowing constraint too loose for unemployed and no edu')
end

if rloop>10
    P            = (1-relax)*(ipc_anddnm*AggY) + relax*Pold;               % Update education price
    Pold         = P;
    P0           = (1-subsidy)*P;                                          % Edu price individuals face
    i_buy        = min(find(gridd.b+abs(gridd.b(1))>=P0));                 % P0 costs i_buy grid points
    i_buymesh2   = gridd.b+abs(gridd.b(1))>=P0;
    x            = gridd.b-P0;
    buildOmega;
end

% Education costs and income tax
flowsINEDUn1 = param.lx1*eU + param.lx2*eEmp + param.lx3_vec(1)*stnHab + param.lx3_vec(2)*stnLab;
flowsINEDUn2 = (eU + eEmp + stn)*param.kappa;
edinflow_old = flowsINEDUn1+flowsINEDUn2;
edcost       = edinflow_old*P;
subs         = subsidy*edcost;
unbenef      = reprate*(unU*gridd.wageL + eU*gridd.wageH); 
switch GT
    case 1
        if rloop>2
            itaxnew(1,1) = unbenef/(gridd.wageL*unEmp + gridd.wageH*eEmp);
            itaxnew(2,1) = itaxnew(1,1) + (subs/(gridd.wageH*eEmp))*(eEmp>=eEmpeps);
            itax         = (1-relax)*itaxnew + relax*itax;
        else
            itax(1,1)    = unbenef/(gridd.wageL*unEmp + gridd.wageH*eEmp);
            itax(2,1)    = itax(1,1) + (subs/(gridd.wageH*eEmp))*(eEmp>=eEmpeps);
        end
        itax        = min(0.5,itax);
        param.itax  = [0 itax(1,1) 0 0 itax(2,1)];
    case 0
        if rloop>2
            itaxnew = (unbenef + subs)/(gridd.wageL*unEmp + gridd.wageH*eEmp);
            itax    = (1-relax)*itaxnew + relax*itax;
        else
            itax    = (unbenef + subs)/(gridd.wageL*unEmp + gridd.wageH*eEmp);
        end
        itax        = min(0.5,itax);
        param.itax  = [0 itax 0 0 itax];
end
for jjj=1:NzNed
    param.itax2(:,jjj) = kron(param.itax(1,jjj),ones(gridd.I,J));
end

% Track evolution of results
MKT              = AggY - AggC - d*KS(rloop) - edcost;
P_r(rloop)       = P;
wr_r(rloop,:)    = gridd.wageH/gridd.wageL;
Y_r(rloop,1)     = AggY;
pop_r(rloop,:)   = [unU,unEmp,stn,eU,eEmp];
tax_r(rloop,:)   = [param.itax(2),param.itax(end)];
MKT_r(rloop,1)   = MKT;
C_r(rloop,:)     = AggC;

% If the free boundary keeps cycling between neighbouring points, pick one (rarely used - increase grid resolution and/or perturb initial guess of Vall to avoid this fix)
switch smallgrid
    case 1
        if rloop>180
           signmovavg = sign(movmean(SS,50));
           countpos   = sum(signmovavg(end-50:end)>0);
           countneg   = sum(signmovavg(end-50:end)<0);
           if (countpos>20)&&(countneg>20)&&(countsmth==0)
               comb      = rloop;
               countsmth = 1;
           end
        end
otherwise
end

if abs(Kdiff(rloop))<(tol_r)&&(rloop>1)&&(abs(MKT)<tol_r)&&(abs(SS(rloop))<tol_r)
%     display('Equilibrium found')
    break
end

end

inequalitystats_SSTSGT;
Res.Welfare = gridd.bdelta'*sum(sum(Vall.*g,3),2);
Res.K       = KD(end);
Res.Y       = AggY;
Res.r       = guess.r;
Res.mktclr  = MKT;
Res.convint = d2;
Res.SSconv  = SS(end);
Res.GDP     = AggY;
Res.tax     = max(itax(:));
Res.rloop   = rloop;
Res.oloop   = oloop;
Res.Meanc   = AggC;
Res.Varc    = Varc;
Res.Giniw   = Giniw;
Res.Ginic   = Ginic;
Res.CVc     = CVc;
Res.CVw     = CVw;
Res.debtors = debtors;
Res.timec   = toc(onerun);
Res.pop     = [unU,unEmp,stn,eU,eEmp];
Res.wage    = [gridd.wageL gridd.wageH];
Res.mnlinc  = Res.pop*gridd.income';
Res.netearn = Res.pop*(gridd.income.*(1-itax))';
Res.b       = gridd.b;
Res.dens    = g;
Res.VnoED   = V;
Res.c       = cons;
Res.UItax   = unbenef/(gridd.wageL*unEmp + gridd.wageH*eEmp);
Res.Vfunc   = Vall;
Res.date    = clock;
Res.P       = P;
Res.aC      = [xi,1-xi];
Res.L       = z_ave;
Res.benchmark_NICLparams = fniniguess;

% Create 'SSresults' folder if it doesn't exist
foldername = fullfile(cd,'SSresults');
if ~exist(foldername,'dir')
    mkdir(foldername);
end

currentDate = datestr(now,'mmmddyyyy'); 
prefname    = fullfile(foldername,currentDate);
filepath    = strcat(prefname, '_bmin_', num2str(abs(round(gridd.b(1), 4) * 10000)),'_I_',num2str(gridd.I),'_');
if GT==0
    matfilenm = strcat(filepath,'TS',num2str(round(100*subsidy)),'GE_','_Edcostscale',num2str(100*round(Edpricescale,2)));
else
    matfilenm = strcat(filepath,'GT',num2str(round(100*subsidy)),'GE_','_Edcostscale',num2str(100*round(Edpricescale,2)));
end
save(matfilenm)