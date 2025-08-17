function vsavehuggett = V_educated_slj(param,guess,num,gridd,Bswitch,vne,lext)

lparams_Vedu; 
% initial guess
if isfield(guess,'vprev')==0
    v0(:,:,1:3) = (1-crra)^(-1)*((gridd.iiincome(:,:,Ned:end).*(1-param.itax2(:,:,Ned:end))...
        + (guess.r+param.kappa*(gridd.bbb(:,:,Ned:end)<0)).*gridd.bbb(:,:,Ned:end)...
        - param.atax(:,:,Ned:end)).^(1-crra)-constant)/(param.rho+param.kappa);                                                                 
else
    v0 = guess.vprev;
end

VstarAux        = zeros(I,J,gridd.Ned);
v               = v0;  

for n=1:maxit 
    V = v;
    % Forward difference for the derivative of v wrt b
    Vbf(1:I-1,:,:)= (V(2:I,:,:)-V(1:I-1,:,:))./gridd.dbbayF(:,:,1:3);
    Vbf(I,:,:)    = (gridd.iiincome(I,:,Ned:end).*(1-param.itax2(I,:,Ned:end)) + guess.r*gridd.bbb(end,:,Ned:end) - param.atax(I,:,Ned:end)).^(-crra); % Boundary condition                                  
    % Backward difference
    Vbb(2:I,:,:)  = (V(2:I,:,:)-V(1:I-1,:,:))./gridd.dbbayB(:,:,1:3);
    Vbb(1,:,:)    = (gridd.iiincome(1,:,Ned:end).*(1-param.itax2(1,:,Ned:end)) + (guess.r+param.kappa*(gridd.bbb(1,:,Ned:end)<0)).*gridd.bbb(1,:,Ned:end) - param.atax(1,:,Ned:end)).^(-crra); % Boundary condition

    Vbf = max(Vbf,10^(-6));                                                % Added for numerical stability
    Vbb = max(Vbb,10^(-6));

    % Consumption and savings with forward difference
    cf = Vbf.^(-1/param.s);  
    sf = gridd.iiincome(:,:,Ned:end).*(1-param.itax2(:,:,Ned:end)) + (guess.r+param.kappa*(gridd.bbb(:,:,Ned:end)<0)).*gridd.bbb(:,:,Ned:end) - cf - param.atax(:,:,Ned:end);
    Hf = (1-param.s)^(-1)*(cf.^(1-param.s)-param.constant) + Vbf.*sf;
    % Consumption and savings with backward difference
    cb = Vbb.^(-1/param.s);
    sb = gridd.iiincome(:,:,Ned:end).*(1-param.itax2(:,:,Ned:end)) + (guess.r+param.kappa*(gridd.bbb(:,:,Ned:end)<0)).*gridd.bbb(:,:,Ned:end) - cb - param.atax(:,:,Ned:end);
    Hb = (1-param.s)^(-1)*(cb.^(1-param.s)-param.constant) + Vbb.*sb;
    % Consumption and derivative of value function at steady state
    c0 = gridd.iiincome(:,:,Ned:end).*(1-param.itax2(:,:,Ned:end)) + (guess.r+param.kappa*(gridd.bbb(:,:,Ned:end)<0)).*gridd.bbb(:,:,Ned:end) - param.atax(:,:,Ned:end);
    H0 = (1-param.s)^(-1)*(c0.^(1-param.s)-param.constant);
    
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

    c        = cf.*If + cb.*Ib + c0.*I0;
    u        = (1-param.s)^(-1)*(c.^(1-param.s)-param.constant);
    % Build A Matrix
    X        = -Ib.*sb./repmat(gridd.dbgridb,1,J,gridd.Ned); % lower
    Y        = -If.*sf./repmat(gridd.dbgridf,1,J,gridd.Ned) + Ib.*sb./repmat(gridd.dbgridb,1,J,gridd.Ned); % center
    Z        = If.*sf./repmat(gridd.dbgridf,1,J,gridd.Ned); % upper
    updiag   = zeros(M,Ned);
    lowdiag  = zeros(M,Ned);
    centdiag = zeros(M,Ned);

    for i = 1:Ned
        centdiag(:,i) = reshape(Y(:,:,i),M,1);
    end

    lowdiag(1:I-1,:) = X(2:I,1,:); 
    updiag(2:I,:) = Z(1:I-1,1,:);  
    for j = 2:J
        lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Ned)]; 
                                                                                   
        updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Ned);squeeze(Z(1:I-1,j,:))];
    end

    AA1               = spdiags(centdiag(:,1),0,M,M)+spdiags(updiag(:,1),1,M,M)+spdiags(lowdiag(:,1),-1,M,M); 
    AA2               = spdiags(centdiag(:,2),0,M,M)+spdiags(updiag(:,2),1,M,M)+spdiags(lowdiag(:,2),-1,M,M);
    AA3               = spdiags(centdiag(:,3),0,M,M)+spdiags(updiag(:,3),1,M,M)+spdiags(lowdiag(:,3),-1,M,M);
    AA                = [AA1, sparse(M,2*M); sparse(M,M), AA2, sparse(M,M); sparse(M,2*M), AA3];
    A                 = AA + Bswitch;                                                  % Adds the transition of other income types
    B                 = spdiags(1/Delta + rho + param.kappa + lext,0,MMe,MMe) - A;
    vec               = u(:) + V(:)/Delta + lext.*vne(:,1) + param.altruism*param.kappa*vne(:,2); 
  
    G1                = griddedInterpolant(bb,aa,V(:,:,1),'makima');
    G2                = griddedInterpolant(bb,aa,V(:,:,2),'makima');
    G3                = griddedInterpolant(bb,aa,V(:,:,3),'makima');
    VstarAux(:,:,1)   = G1(gridd.btilde,gridd.atilde);
    VstarAux(:,:,2)   = G2(gridd.btilde,gridd.atilde);
    VstarAux(:,:,3)   = G3(gridd.btilde,gridd.atilde);  
    VstarAux(:,end,:) = -5000;
    VstarAux(1,:,:)   = -5000;
    VstarAux_stacked  = VstarAux(:);
    q                 = -vec + B*VstarAux_stacked;
    z0                = V(:)-VstarAux_stacked; 
    zvec              = LCP(B,q,zeros(gridd.MMe,1),Inf*ones(gridd.MMe,1),z0,0); % Solve for V
    LCP_error         = max(abs(zvec.*(B*zvec + q)));
    V_stacked_123     = zvec+VstarAux_stacked;
    V                 = reshape(V_stacked_123,I,J,3);
    Vchange           = V - v;
    v                 = V;
    dist(n)           = max(max(max(abs(Vchange))));
    if dist(n)<tol_hjb
        break
    end
end

adothuggett             = sf.*If + sb.*Ib; 
vsavehuggett.V          = V;
vsavehuggett.c          = c;
vsavehuggett.adot       = adothuggett;
vsavehuggett.Amat       = A;
vsavehuggett.RepayEarly = abs(V_stacked_123 - VstarAux_stacked)<param.etol;
vsavehuggett.distn      = dist;
end