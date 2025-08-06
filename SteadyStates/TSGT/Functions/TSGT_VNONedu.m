function vsavehuggett = TSGT_VNONedu(param,guess,num,gridd,Aswitch,bnodebt)

translateparam_nonedu;
% Initial guess
v0(:,1:Nz) = (1-param.s)^(-1)*((gridd.iiincome(:,1:Nz).*(1-param.itax2(:,1:Nz)) +...
             (guess.r+param.kappa*(gridd.bb(:,1:Nz)<0)).*gridd.bb(:,1:Nz)).^(1-param.s)-param.constant)/(param.rho+param.kappa);  
v = v0;    
for n=1:num.maxit_hjb
    V = v;
    % Forward difference for the derivative of v wrt b
    Vbf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))./gridd.dbyF(:,1:Nz); 
    Vbf(I,:)     = (gridd.iiincome(I,1:Nz).*(1-param.itax2(I,1:Nz)) + guess.r*gridd.bmax ).^(-param.s); % Boundary condition
    % Backward difference
    Vbb(2:I,:)   = (V(2:I,:)-V(1:I-1,:))./gridd.dbyB(:,1:Nz);
    Vbb(1,:)     = (gridd.iiincome(1,1:Nz).*(1-param.itax2(1,1:Nz)) + (guess.r+param.kappa*(gridd.b(1)<0))*gridd.b(1)).^(-param.s); % Boundary condition

    Vbf = max(Vbf,10^(-6));                                                % For numerical stability  
    Vbb = max(Vbb,10^(-6));

    % Consumption and savings with forward difference
    cf         = Vbf.^(-1/param.s); 
    sf         = gridd.iiincome(:,1:Nz).*(1-param.itax2(:,1:Nz)) + (guess.r+param.kappa*(gridd.bb(:,1:Nz)<0)).*gridd.bb(:,1:Nz) - cf;
    Hf         = (1-param.s)^(-1)*(cf.^(1-param.s)-param.constant) + Vbf.*sf;
    % Consumption and savings with backward difference
    cb         = Vbb.^(-1/param.s); 
    sb         = gridd.iiincome(:,1:Nz).*(1-param.itax2(:,1:Nz)) + (guess.r+param.kappa*(gridd.bb(:,1:Nz)<0)).*gridd.bb(:,1:Nz) - cb;
    Hb         = (1-param.s)^(-1)*(cb.^(1-param.s)-param.constant) + Vbb.*sb;
    % Consumption and derivative of value function at steady state
    c0         = gridd.iiincome(:,1:Nz).*(1-param.itax2(:,1:Nz)) + (guess.r+param.kappa*(gridd.bb(:,1:Nz)<0)).*gridd.bb(:,1:Nz);
    H0         = (1-param.s)^(-1)*(c0.^(1-param.s)-param.constant);
   
    if param.AssumeConcavity == 1
        Ib = (sb < 0); 
        If = (sf > 0).*(1-Ib); 
        I0 = (1-If-Ib);
    else 
        Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0);
        Iboth = (sb<0).*(sf>0);
        Ib = Iunique.*(sb<0).*(Hb>H0) + Iboth.*(Hb>Hf).*(Hb>H0);
        If = Iunique.*(sf>0).*(Hf>H0) + Iboth.*(Hf>Hb).*(Hf>H0);
        I0 = 1-Ib-If;
    end

    c         = cf.*If + cb.*Ib + c0.*I0;
    u         = (1-param.s)^(-1)*(c.^(1-param.s)-param.constant);
    X         = -Ib.*sb./repmat(gridd.dbgridb,1,gridd.Nz); % lower
    Y         = -If.*sf./repmat(gridd.dbgridf,1,gridd.Nz) + Ib.*sb./repmat(gridd.dbgridb,1,gridd.Nz); % center
    Z         = If.*sf./repmat(gridd.dbgridf,1,gridd.Nz);  % upper
    A1        = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2        = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A         = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch;
    B         = (1/num.Delta + param.rho + param.kappa)*speye(I*Nz) - A;
    u_stacked = reshape(u,I*gridd.Nz,1);
    V_stacked = reshape(V,I*gridd.Nz,1);
    vec       = u_stacked + V_stacked/num.Delta + param.altruism*param.kappa*repmat(bnodebt*V(:,1),2,1);
    V_stacked = B\vec;
    V         = reshape(V_stacked,I,gridd.Nz);
    Vchange   = V - v;
    v         = V;
    dist(n)   = max(max(abs(Vchange)));
    if dist(n)<num.tol_hjb
%         sprintf('Value function converged')
        break
    end
end
adothuggett       = sf.*If + sb.*Ib; 
vsavehuggett.V    = V;
vsavehuggett.c    = c;
vsavehuggett.adot = adothuggett;
vsavehuggett.Amat = A;
vsavehuggett.dist = dist(n);
vsavehuggett.nloop= n;
end