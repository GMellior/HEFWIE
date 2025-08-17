
for sss=1:abtypes
    blocklx     = [[spdiags(param.lx3_vec(sss)*ones(M,1),0,M,M),sparse(M,M)];spdiags([param.lx1*ones(M,1);param.lx2*ones(M,1)],0,MMu,MMu)];
    adjpoints   = pointslist(logical(adjMMM2(:,sss)+adjMMM2_ER(:,sss)));
    noadjpoints = pointslist(~logical(adjMMM2(:,sss)+adjMMM2_ER(:,sss)));
    % adj already has i_buymesh2 -  and also should only pick adjMMM2 and adjMMM2_ER
    bAdjAux     = repmat(bAdj1(:),NzNed,1).*adjMMM2(:,sss) + adjMMM2_ER(:,sss).*repmat(gridd.btilde(:),NzNed,1); 
    aAdjAux     = repmat(aAdj1(:),NzNed,1).*adjMMM2(:,sss) + adjMMM2_ER(:,sss).*repmat(gridd.atilde(:),NzNed,1); 
    shift       = 0;
    % BUILD M MATRIX - put 1's on diagonal for points not in adjustment region
    MM          = sparse(noadjpoints, noadjpoints, ones(length(noadjpoints),1), MMM, MMM);
for mcount = 1:length(adjpoints)
    m            = adjpoints(mcount);
    bprime       = bAdjAux(m,1);
    aprime       = aAdjAux(m,1);
    if m <= M 
        mm       = (Nz+1)*adjMMM2(m,sss) + 1*adjMMM2_ER(m,sss);
    elseif (m>M)&&(m<=MMu)
        mm       = (Nz+1)*adjMMM2(m,sss) + 2*adjMMM2_ER(m,sss);
    else
        mm       = ceil(m/M);
    end
    bprime_left  = discretize(bprime, gridd.b); % Index of the left edge for cash
    aprime_left  = discretize(aprime, gridd.a); % Student loans
    bprime_right = bprime_left + 1;
    aprime_right = aprime_left + 1;
    point11      = pointsmat(bprime_left, aprime_left,  mm);
    point12      = pointsmat(bprime_left, aprime_right, mm);
    point21      = pointsmat(bprime_right, aprime_left, mm);
    point22      = pointsmat(bprime_right, aprime_right,mm);
    neighpoints  = [point11 point12;
                    point21 point22];
    totarea = (gridd.b(bprime_right) - gridd.b(bprime_left)) * (gridd.a(aprime_right) - gridd.a(aprime_left)); % this is just da * db with equispaced grid 
    weights = totarea^(-1) * [gridd.b(bprime_right) - bprime; bprime - gridd.b(bprime_left)] * [gridd.a(aprime_right) - aprime aprime - gridd.a(aprime_left)];
    MM      = MM + sparse(m*ones(4, 1), shift+neighpoints(:), weights(:), MMM, MMM);
end

% Eliminate falling into adjustment points
for j=1:100
    MMn  = MM*MM;
    diff = max(abs(MMn-MM));
    MM   = MMn;
    if diff<1e-9
        % disp('MM Converged')
        break
    end
end

MMsss{sss} = MM;
Am{sss}    = [ANC{sss},sparse(MMu,MMe);blocklx,AC{sss}-spdiags([param.lx3_vec(sss)*ones(M,1);param.lx1*ones(M,1);param.lx2*ones(M,1)],0,MMe,MMe)]...
              -spdiags(param.kappa*ones(MMM,1),0,MMM,MMM) + param.kappa*prob_stay(:,sss).*born_rentry;

end

Afinal = [[Am{1}*MMsss{1},(param.kappa*(1-prob_stay(:,1)).*spdiags(1*~logical(adjMMM2(:,1)+adjMMM2_ER(:,1)),0,MMM,MMM)*born_rentry)*MMsss{2}];...
         [(param.kappa*(1-prob_stay(:,2)).*spdiags(1*~logical(adjMMM2(:,2)+adjMMM2_ER(:,2)),0,MMM,MMM)*born_rentry)*MMsss{1},Am{2}*MMsss{2}]];

adjpoints1               = pointslist(logical(adjMMM2(:,1)+adjMMM2_ER(:,1)));
adjpoints2               = pointslist(logical(adjMMM2(:,2)+adjMMM2_ER(:,2)));
Afinal(adjpoints1,:)     = 0;
Afinal(MMM+adjpoints2,:) = 0;
Afinal                   = Afinal - constantkfe2*[[spdiags(adjMMM2(:,1)+adjMMM2_ER(:,1),0,MMM,MMM),sparse(MMM,MMM)];[sparse(MMM,MMM),spdiags(adjMMM2(:,2)+adjMMM2_ER(:,2),0,MMM,MMM)]];
Afinal                   = Afinal + constantkfe2*[[(MMsss{1}-spdiags(1*~logical(adjMMM2(:,1)+adjMMM2_ER(:,1)),0,MMM,MMM)),sparse(MMM,MMM)];[sparse(MMM,MMM),(MMsss{2}-spdiags(1*~logical(adjMMM2(:,2)+adjMMM2_ER(:,2)),0,MMM,MMM))]];
AfT                      = Afinal';
b                        = [0.1;zeros(MMMM-1,1)];
AfT(1,:)                 = [0.1,zeros(1,MMMM-1)];
g_tilde                  = AfT\b;
g_sum                    = g_tilde'*ones(MMMM,1);
g_tilde                  = g_tilde./g_sum;
gg                       = grid_diag\g_tilde; % Convert from \tilde{g} to g
g                        = reshape(gg,gridd.I,J,NzNed,abtypes);