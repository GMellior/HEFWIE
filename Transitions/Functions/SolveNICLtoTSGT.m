
% Prices and tax rate pathsloop
for rloop=1:2000
    Vall            = v_st;
    % Time iteration of the value function
for t=perds:-1:1
    update_tran_prices_TSGT;
    Omega_t{t} = Omega;
    for sss=1:abtypes
        if t==perds
            guess.vprev          = v_st(:,Nz+1:end,sss);
        else
            guess.vprev          = Vall(:,Nz+1:end,sss);
        end
        param.lx3                = param.lx3_vec(sss);
        computeVs_SFTSGT_2025_transition;
        cons(:,1:Nz,sss)         = c;
        cons(:,Nz+1:NzNed,sss)   = v0pack2.c;
        adj                      = (abs(V_stacked - Vstar_stacked)<etol);   % Allow for some tolerance in defining solution for adjustment decision
        adjmat                   = reshape(adj,I,Nz);
        adjmat2                  = zeros(I,Nz);
        adjmat0(:,1)             = (adjmat(:,1)*1).*indxpts;
        adjmat0(:,2)             = (adjmat(:,2)*1).*indxpts;
        adjmat0(adjmat0==0)      = I+1;
        for jjj=1:J
            [row(jjj,:),~]    = min(adjmat0);
            [rowmax(jjj,:),~] = max(adjmat0);
        end
        if (rloop>50)&&(rloop<1800)                                        % See comments on steady state file
            for jjj=1:J
                adj           = round( (1/(comb3+1))*(adj.*repmat(i_buymesh2(:),2,1) + squeeze(sum(adj_tr(:,end,sss,t),2)) ));
            end
        elseif rloop>1799
            adj               = round( (1/(comb4+1))*(adj.*repmat(i_buymesh2(:),2,1) + squeeze(sum(adj_tr(:,1:comb4,sss,t),2)) ));
        end
        adj                      = logical(adj.*repmat(i_buymesh2(:),2,1)); % Added to avoid choosing edu at forbidden area
        adjMMM                   = [adj;false(MMe,1)];
        adjRegion1               = reshape(adj(1:M),I,J).*ones(I,J);        % Multiplication by ones(I,J) transforms logicals into floats
        adjRegion2               = reshape(adj(M+1:MMu),I,J).*ones(I,J);
        adjRegion(:,1)           = adjRegion1; 
        adjRegion(:,2)           = adjRegion2; 
        adjMMM2(:,sss)           = adjMMM;
        adj_tr(:,1:end-1,sss,t)  = adj_tr(:,2:end,sss,t);
        adj_tr(:,end,sss,t)      = adj;
        ANC{sss}                 = A;
        AC{sss}                  = v0pack2.Amat;
        Vstar_stacked_sss(:,sss) = Vstar_stacked;
    end
    %%
    getdistribution_TSGT_transition;
    V_t(:,:,:,t)    = Vall;
    cons_t(:,:,:,t) = cons;
end
%%%%% KFE
gold        = gg0;
for n=1:perds
    AT        = C_t{n}';
    num.Delta = dtsvec(n);
    compute_edcost_Amatrix_tran;
    if n==1
        gg{n}   = (speye(MMMM) - AT*num.Delta)\(gg0.*bayyydelta(:));
        g       = full(reshape(ndSparse(grid_diag\gg{n}),I,NzNed,abtypes));
    else
        gg{n}   = (speye(MMMM) - AT*num.Delta)\gg{n-1};
        g       = full(reshape(ndSparse(grid_diag\gg{n}),I,NzNed,abtypes));
    end

    gold        = g(:);
    K_t(n,1)    = (gridd.b.*gridd.bdelta)'*sum(sum(g,3),2) + pathA(n);
    %%% Population groups
    unU_t(n,1)    = gridd.bdelta'*sum(g(:,1,:),3);                         % Share unemployed with no education
    unEmp_t(n,1)  = gridd.bdelta'*sum(g(:,2,:),3);                         % Share employed with no education
    stn_t(n,1)    = gridd.bdelta'*sum(g(:,3,:),3);                         % Share student
    eU_t(n,1)     = gridd.bdelta'*sum(g(:,4,:),3);                         % Share unemployed WITH education
    eEmp_t(n,1)   = gridd.bdelta'*sum(g(:,5,:),3);                         % Share employed WITH education
    stnHab_t(n,1) = gridd.bdelta'*g(:,3,1);                                % Share student high ability (low dropout)
    stnLab_t(n,1) = gridd.bdelta'*g(:,3,2);                                % Share student low ability
    AggC_t(n,1)   = gridd.bdelta'*(sum(sum(g.*cons_t(:,:,:,n),3),2));      % Aggregate consumption
    inequalitystats_TSGT_transition;                                       % Gini and debtor shares
end
% Update price and tax rates paths
z_ave_t            = (xi0*((unEmp_t+zstn*stn_t).^nu) + (1-xi0)*((max(eEmp_t,eEmpeps)).^nu)).^(1/nu);
rnew_t             = ((z_ave_t).^(1-param.alpha)).*(param.alpha* Aprod*K_t.^(param.alpha-1)) - d;          % Interest rates
wL_tn              = min(((1-param.alpha)*Aprod)*((K_t.^param.alpha).*(z_ave_t.^(1-param.alpha-nu))).*(xi0*(unEmp_t+zstn*stn_t).^(nu-1)),1);
wH_tn              = min(max(((1-param.alpha)*Aprod)*((K_t.^param.alpha).*(z_ave_t.^(1-param.alpha-nu))).*((1-xi0)*(max(eEmp_t,eEmpeps)).^(nu-1)),wL_tn),2*wL_tn);
wL_n_t             = wL_tn;
wH_n_t             = wH_tn;
r_t                = relax.*r_t  + (1-relax).*rnew_t;
wL_t               = relax.*wL_t + (1-relax).*wL_n_t;
wH_t               = relax.*wH_t + (1-relax).*wH_n_t;
Knew_t(:,rloop+1)  = K_t;                  
checkconv          = max(abs(K_t-Knew_t(:,rloop)));
checkconv_r(rloop) = checkconv;            %#ok

if rloop>95
    P_t         = relax.*P_told + (1-relax).*(ipc_anddnm*Aprod*(K_t.^param.alpha).*(z_ave_t.^(1-param.alpha)));
    P_told      = P_t;
end

if mod(rloop,1)==0 % Increse grid resolution or smooth bumps in tax rate - see steady state file comments
    edinflow_t_old(160:end)    = 0.4*hpfilter(edinflow_t_old(160:end),100) + 0.6*edinflowT;
    edinflow_t(160:end)        = 0.4*hpfilter(edinflow_t(160:end),100) + 0.6*edinflowT;
    edinflow_t_old(end-15:end) = 0.05*edinflow_t_old(end-15:end) + 0.95*edinflowT;
    edinflow_t(end-15:end)     = 0.05*edinflow_t(end-15:end) + 0.95*edinflowT;
    eo                         = csaps(Time,edinflow_t,0.05);
    edinflow_t                 = fnval(eo,Time);
end

edinflow_t(end)     = edinflow;
edinflow_t_old      = edinflow_t;     
edinflow_t_old(end) = edinflow;
edc1                = subsidy*edinflow_t_old.*P_t;
unbenef             = reprate*(unU_t.*wL_t + eU_t.*wH_t);
itax_t0             = (unbenef - pathA.*(amort+r_t))./(wL_t.*unEmp_t + wH_t.*eEmp_t);

switch GT
    case 0
        itax        = itax_t0 + edc1./(wL_t.*unEmp_t + wH_t.*eEmp_t);
        itax_t      = min(0.9,itax);
        itax_t      = relax.*itax_told + (1-relax).*itax_t;
    case 1
        itax_t      = [itax_t0,...
                       itax_t0 + (edc1./(wH_t.*eEmp_t)).*(eEmp_t>=eEmpeps)];
        itax_t(:,1) = relax.*itax_told(:,1) + (1-relax).*itax_t(:,1);
        itax_t(:,2) = relax.*itax_told(:,2) + (1-relax).*itax_t(:,2);
end

itax_told = itax_t;

switch GT
    case 0
        meaninc_t = wL_t.*(reprate*unU_t + (1-itax_t).*unEmp_t + zstn*stn_t) + wH_t.*(reprate*eU_t + (1-itax_t).*eEmp_t); % using old wage     
    case 1
        meaninc_t = wL_t.*(reprate*unU_t + (1-itax_t).*unEmp_t + zstn*stn_t) + wH_t.*(reprate*eU_t + (1-itax_t).*eEmp_t); % using old wage     
end

AggY_t    = Aprod*(K_t.^param.alpha).*(z_ave_t.^(1-param.alpha));
AggKdot_t = AggY_t - AggC_t - d*K_t - edinflow_t_old.*P_t;

if checkconv<1e-6
    break
end

end