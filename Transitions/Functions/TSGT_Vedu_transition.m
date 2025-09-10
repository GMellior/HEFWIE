function vres = TSGT_Vedu_transition(param,guess,num,gridd,Aswitch,vne,lext)
translateparam_edu;
% Initial guess
if isfield(guess,'vprev')==0
    v(:,1:Ned) = (1-param.s)^(-1)*((gridd.iiincome(:,Nz+1:end).*(1-param.itax2(:,Nz+1:end)) +...
                   (guess.r+param.kappa*(gridd.bb(:,Nz+1:end)<0)).*gridd.bb(:,Nz+1:end)).^(1-param.s)-param.constant)/(param.rho+param.kappa); 
else
    v = guess.vprev;
end

V_stacked = v(:);
for n=1:10
    V            = v;
    Vbf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))./gridd.dbyF;                     % Forward difference for the derivative of v wrt b
    Vbf(I,:)     = (gridd.iiincome(I,Nz+1:end).*(1-param.itax2(I,Nz+1:end)) + guess.r*gridd.bmax ).^(-param.s);
    Vbb(2:I,:)   = (V(2:I,:)-V(1:I-1,:))./gridd.dbyB;                     % Backward difference
    Vbb(1,:)     = (gridd.iiincome(1,Nz+1:end).*(1-param.itax2(1,Nz+1:end)) + (guess.r+param.kappa*(gridd.b(1)<0))*gridd.b(1)).^(-param.s); 

    Vbf          = max(Vbf,10^(-6)); 
    Vbb          = max(Vbb,10^(-6));
    cf           = Vbf.^(-1/param.s);                                       % Consumption and savings with forward difference
    sf           = gridd.iiincome(:,Nz+1:end).*(1-param.itax2(:,Nz+1:end)) + (guess.r+param.kappa*(gridd.bb(:,Nz+1:end)<0)).*gridd.bb(:,Nz+1:end) - cf;
    Hf           = (1-param.s)^(-1)*(cf.^(1-param.s)-param.constant) + Vbf.*sf;
    cb           = Vbb.^(-1/param.s);                                       % Consumption and savings with backward difference
    sb           = gridd.iiincome(:,Nz+1:end).*(1-param.itax2(:,Nz+1:end)) + (guess.r+param.kappa*(gridd.bb(:,Nz+1:end)<0)).*gridd.bb(:,Nz+1:end) - cb;
    Hb           = (1-param.s)^(-1)*(cb.^(1-param.s)-param.constant) + Vbb.*sb;
    c0           = gridd.iiincome(:,Nz+1:end).*(1-param.itax2(:,Nz+1:end)) + (guess.r+param.kappa*(gridd.bb(:,Nz+1:end)<0)).*gridd.bb(:,Nz+1:end);
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

    c = cf.*If + cb.*Ib + c0.*I0;
    u = (1-param.s)^(-1)*(c.^(1-param.s)-param.constant);

    % Build A matrix
    X         = -Ib.*sb./repmat(gridd.dbgridb,1,gridd.Ned);                % lower
    Y         = -If.*sf./repmat(gridd.dbgridf,1,gridd.Ned) + Ib.*sb./repmat(gridd.dbgridb,1,gridd.Ned); % center
    Z         = If.*sf./repmat(gridd.dbgridf,1,gridd.Ned);                 % upper
    A1        = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
    A2        = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
    A3        = spdiags(Y(:,3),0,I,I)+spdiags(X(2:I,3),-1,I,I)+spdiags([0;Z(1:I-1,3)],1,I,I);
    A         = [A1,sparse(I,2*I);sparse(I,I),A2,sparse(I,I);sparse(I,2*I),A3] + Aswitch;
    B         = (1/num.Delta + param.rho + lext + param.kappa).*speye(I*Ned) - A;
    u_stacked = reshape(u,I*gridd.Ned,1);
    vec       = u_stacked + V_stacked/num.Delta + lext.*vne(:,1) + param.altruism*param.kappa*vne(:,2);
    % New V
    V_stacked_123 = B\vec;
    V             = reshape(V_stacked_123,I,gridd.Ned);
    Vchange       = V - v;
    v             = V;
    dist(n)       = max(max(abs(Vchange)));
    if dist(n)<num.tol_hjb
        break
    end

end
vres.V    = V;
vres.c    = c;
vres.adot = sf.*If + sb.*Ib;
vres.Amat = A;
vres.dist = dist(n);
vres.nloop= n;
end