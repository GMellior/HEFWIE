
for sss=1:abtypes
    blocklx     = [[spdiags(param.lx3_vec(sss)*ones(M,1),0,M,M),sparse(M,M)];spdiags([param.lx1*ones(M,1);param.lx2*ones(M,1)],0,MMu,MMu)];
    adjpoints   = pointslist(logical(adjMMM2(:,sss)));
    noadjpoints = pointslist(~logical(adjMMM2(:,sss)));
    % adj already has i_buymesh2 -  and also should only pick adjMMM2 and adjMMM2_ER
    bAdjAux     = repmat(gridd.b-(1-subsidy)*P_t(t),NzNed,1).*adjMMM2(:,sss); 
    % BUILD M MATRIX - put 1's on diagonal for points not in adjustment region
    MM          = sparse(noadjpoints, noadjpoints, ones(length(noadjpoints),1), MMM, MMM);
for mcount = 1:length(adjpoints)
    m            = adjpoints(mcount);
    bprime       = bAdjAux(m,1);
    if m <= M 
        mm       = (Nz+1)*adjMMM2(m,sss);
    elseif (m>M)&&(m<=MMu)
        mm       = (Nz+1)*adjMMM2(m,sss);
    else
        mm       = ceil(m/M);
    end
    bprime_left  = discretize(bprime, gridd.b);
    bprime_right = bprime_left + 1;
    point11      = pointsmat(bprime_left,mm);
    point12      = pointsmat(bprime_right,mm);
    neighpoints  = [point11 point12];
    totarea      = (gridd.b(bprime_right) - gridd.b(bprime_left));
    weights      = totarea^(-1) * [gridd.b(bprime_right) - bprime bprime - gridd.b(bprime_left)];
    MM           = MM + sparse(m*ones(2, 1), neighpoints(:), weights(:), MMM, MMM);
end

% Eliminate falling into adjustment points
for j=1:100
    MMn  = MM*MM;
    diff = max(abs(MMn-MM));
    MM   = MMn;
    if diff<1e-9
        break
    end
end

MMsss{sss} = MM;
Am{sss}    = [ANC{sss},sparse(MMu,MMe);blocklx,AC{sss}-spdiags([param.lx3_vec(sss)*ones(M,1);param.lx1*ones(M,1);param.lx2*ones(M,1)],0,MMe,MMe)]...
              -spdiags(param.kappa*ones(MMM,1),0,MMM,MMM) + param.kappa*prob_stay(:,sss).*born_rentry;

end

Afinal                   = [[Am{1}*MMsss{1},(param.kappa*(1-prob_stay(:,1)).*spdiags(1*~logical(adjMMM2(:,1)),0,MMM,MMM)*born_rentry)*MMsss{2}];...
                           [(param.kappa*(1-prob_stay(:,2)).*spdiags(1*~logical(adjMMM2(:,2)),0,MMM,MMM)*born_rentry)*MMsss{1},Am{2}*MMsss{2}]];
adjpoints1               = pointslist(logical(adjMMM2(:,1)));
adjpoints2               = pointslist(logical(adjMMM2(:,2)));
Afinal(adjpoints1,:)     = 0;
Afinal(MMM+adjpoints2,:) = 0;
Afinal                   = Afinal - [[spdiags(adjMMM2(:,1),0,MMM,MMM),sparse(MMM,MMM)];[sparse(MMM,MMM),spdiags(adjMMM2(:,2),0,MMM,MMM)]]*dtsvec(t);
Afinal                   = Afinal + [[(MMsss{1}-spdiags(1*~logical(adjMMM2(:,1)),0,MMM,MMM)),sparse(MMM,MMM)];[sparse(MMM,MMM),(MMsss{2}-spdiags(1*~logical(adjMMM2(:,2)),0,MMM,MMM))]]*dtsvec(t);
C_t{t}                   = Afinal;
laction{t}               = adjMMM2;