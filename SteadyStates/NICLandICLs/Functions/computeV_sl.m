
for oloop=1:20 %% Outer HJB loop
    lext           = kron([param.lx3;param.lx1;param.lx2],ones(M,1));
    vne(:,1)       = [kron(ones(Nz,1),reshape(Vall(:,:,1,sss),M,1));reshape(Vall(:,:,2,sss),M,1)];
    vne(1:M,2)     = reshape(bnodebt*(  (1-parent_abs_small(1,sss))*repmat(squeeze(max(Vall(:,end,1,sss),Vall(gridd.bzeroind+1,end,1,sss))),1,J)+ parent_abs_small(1,sss)*repmat(squeeze(max(Vall(:,end,1,mod(sss,2)+1),Vall(gridd.bzeroind+1,end,1,mod(sss,2)+1))),1,J)),M,1);
    vne(M+1:end,2) = repmat(reshape(bnodebt*(  (1-parent_abs_small(2,sss))*repmat(squeeze(max(Vall(:,end,1,sss),Vall(gridd.bzeroind+1,end,1,sss))),1,J)+ parent_abs_small(2,sss)*repmat(squeeze(max(Vall(:,end,1,mod(sss,2)+1),Vall(gridd.bzeroind+1,end,1,mod(sss,2)+1))),1,J)),M,1),gridd.Nz,1);
    v0pack2        = V_educated_sl(param,guess,num,gridd,Bswitch2,vne,lext); % Inner loop V_edu
    Vall(:,:,Nz+1:end,sss) = v0pack2.V;
    v                      = Vall(:,:,1:Nz,sss);
for n=1:maxit %% Inner loop - solve for value of not educated
    V = v;
    % Forward difference for the derivative of v wrt b
    Vbf(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))./gridd.dbbayF(:,:,1:Nz);
    Vbf(I,:,:)     = (gridd.iiincome(I,:,1:Nz).*(1-param.itax2(I,:,1:Nz)) + guess.r*gridd.bbb(I,:,1:Nz) - param.atax(I,:,1:Nz)).^(-crra); % Boundary condition
    % Backward difference
    Vbb(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))./gridd.dbbayB(:,:,1:Nz);
    Vbb(1,:,:)   = (gridd.iiincome(1,:,1:Nz).*(1-param.itax2(1,:,1:Nz)) + (guess.r+param.kappa*(gridd.bbb(1,:,1:Nz)<0)).*gridd.bbb(1,:,1:Nz) - param.atax(1,:,1:Nz)).^(-crra); % Boundary condition

    Vbf = max(Vbf,10^(-6));                                                % For numerical stability
    Vbb = max(Vbb,10^(-6));

    % Consumption and savings with forward difference               
    cf         = Vbf.^(-1/param.s); 
    sf         = gridd.iiincome(:,:,1:Nz).*(1-param.itax2(:,:,1:Nz)) + (guess.r+param.kappa*(gridd.bbb(:,:,1:Nz)<0)).*gridd.bbb(:,:,1:Nz) - cf - param.atax(:,:,1:Nz);
    Hf         = (1-param.s)^(-1)*(cf.^(1-param.s)-param.constant) + Vbf.*sf;
    % Consumption and savings with backward difference
    cb         = Vbb.^(-1/param.s); 
    sb         = gridd.iiincome(:,:,1:Nz).*(1-param.itax2(:,:,1:Nz)) + (guess.r+param.kappa*(gridd.bbb(:,:,1:Nz)<0)).*gridd.bbb(:,:,1:Nz) - cb - param.atax(:,:,1:Nz);
    Hb         = (1-param.s)^(-1)*(cb.^(1-param.s)-param.constant) + Vbb.*sb;
    % Consumption and derivative of value function at temporary steady state
    c0         = gridd.iiincome(:,:,1:Nz).*(1-param.itax2(:,:,1:Nz)) + (guess.r+param.kappa*(gridd.bbb(:,:,1:Nz)<0)).*gridd.bbb(:,:,1:Nz) - param.atax(:,:,1:Nz);
    indicsubs1 = (c0(:,:,1)<10^(-6));
    sb(1,:,1)  = max(sb(1,:,1),10^(-6));
    H0         = (1-param.s)^(-1)*(c0.^(1-param.s)-param.constant);
   
    if param.AssumeConcavity == 1
        Ib      = (sb < 0); 
        If      = (sf > 0).*(1-Ib); 
        I0      = (1-If-Ib);
    else 
        Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0);
        Iboth   = (sb<0).*(sf>0);
        Ib      = Iunique.*(sb<0).*(Hb>H0) + Iboth.*(Hb>Hf).*(Hb>H0);
        If      = Iunique.*(sf>0).*(Hf>H0) + Iboth.*(Hf>Hb).*(Hf>H0);
        I0      = 1-Ib-If;
    end

    c = cf.*If + cb.*Ib + c0.*I0;
    u = (1-param.s)^(-1)*(c.^(1-param.s)-param.constant);


    X = -Ib.*sb./repmat(gridd.dbgridb,1,J,gridd.Nz); % lower
    Y = -If.*sf./repmat(gridd.dbgridf,1,J,gridd.Nz) + Ib.*sb./repmat(gridd.dbgridb,1,J,gridd.Nz); % center
    Z = If.*sf./repmat(gridd.dbgridf,1,J,gridd.Nz);  % upper

    updiag   = zeros(M,Nz);
    lowdiag  = zeros(M,Nz);
    centdiag = zeros(M,Nz);

    for i = 1:Nz
        centdiag(:,i) = reshape(Y(:,:,i),M,1);
    end
    lowdiag(1:I-1,:) = X(2:I,1,:); 
    updiag(2:I,:)    = Z(1:I-1,1,:); 
    for j = 2:J
        lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];  
        updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
    end
    AA1=spdiags(centdiag(:,1),0,M,M)+spdiags(updiag(:,1),1,M,M)+spdiags(lowdiag(:,1),-1,M,M); 
    AA2        = spdiags(centdiag(:,2),0,M,M)+spdiags(updiag(:,2),1,M,M)+spdiags(lowdiag(:,2),-1,M,M);
    AA         = [AA1, sparse(M,M); sparse(M,M), AA2];
    A          = AA + Bswitch; 
    B          = (1/num.Delta + rho + param.kappa)*speye(MMu) - A;
    u_stacked  = reshape(u,MMu,1);
    V_stacked  = reshape(V,MMu,1);
    Vstudent   = griddedInterpolant(gridd.bb,gridd.aa,v0pack2.V(:,:,1),'spline');
    VEarlyRep1 = griddedInterpolant(gridd.bb,gridd.aa,V(:,:,1));
    VEarlyRep2 = griddedInterpolant(gridd.bb,gridd.aa,V(:,:,2));

%%%%%%% Vstar - option value of being educated but having less assets (because of edu cost)
    % Create gridded interpolant for value of being a student
    VstarAuxvec     = Vstudent(bAdjvec1(:,:,1),aAdjvec1(:,:,1));
    VstarAux        = reshape(VstarAuxvec,I,J,1);
    VstarAux(:,:,1) = (-5000 + VstarAux(:,:,1)).*(sumGridP==-5000) + VstarAux(:,:,1).*(sumGridP~=-5000)...
                      - param.mcost(sss)./(1+0*gridd.bb.*(gridd.bb>=0));
    VstarAux(:,:,2) = VstarAux(:,:,1);
    Vstar_GoToUni  = reshape(VstarAux,MMu,1);
%%%%%%% Vstar - option value of early repayment of student loans
    VstarER(:,:,1)   = VEarlyRep1(gridd.btilde,gridd.atilde);
    VstarER(:,:,2)   = VEarlyRep2(gridd.btilde,gridd.atilde);
    VstarER(:,end,:) = -5000;
    VstarER(1,:,:)   = -5000;
    VstarAll_stacked = max(VstarER(:),Vstar_GoToUni);
    veclcp           = u_stacked + V_stacked/num.Delta  + 0.001*max(VstarAll_stacked-V_stacked,0)...
                       + param.altruism*param.kappa*repmat(reshape(bnodebt*( (1-parent_abs_small(1,sss))*repmat(squeeze(Vall(:,end,1,sss)),1,J)+parent_abs_small(1,sss)*repmat(squeeze(Vall(:,end,1,mod(sss,2)+1)),1,J)),M,1),gridd.Nz,1);
    q                = -veclcp + B*VstarAll_stacked;
    z0               = V_stacked-VstarAll_stacked; 
    zvec             = LCP(B,q,lbnd,ubnd,z0,0);         % Solve for V
    LCP_error        = max(abs(zvec.*(B*zvec + q)));
    V_stacked        = zvec+VstarAll_stacked;
    V(:,:,1)         = reshape(V_stacked(1:M),I,J);
    V(:,:,2)         = reshape(V_stacked(M+1:MMu),I,J);
    Vchange          = V - v;
    v                = V;

    Vall(:,:,1:Nz,sss) = V;
    LCP_errorchange  = LCP_error - LCP_errorold;
    LCP_errorold     = LCP_error;
    if LCP_errorchange==0 % Sometimes the LCP error gets stuck
       LCP_error     = 1e-10; 
    end
    dist(n) = max([max(abs(Vchange(:))) LCP_error/100]); % Make sure LCP error is small
    if (dist(n)<num.tol_hjb)
        break
    end
end
    d1             = dist(n);
    vnew           = v;
    vdistol(oloop) = max(max(max(abs(vnew-vold))));
    vold           = vnew;
    guess.vprev    = v0pack2.V; 
    d2             = vdistol(oloop);
    if vdistol(oloop)<num.tol_hjb
%         disp('Outer V loop converged')
        break
    end
end