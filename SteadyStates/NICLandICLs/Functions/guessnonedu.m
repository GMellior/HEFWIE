
for n=1:maxit
    V = v;
    % Forward difference for the derivative of v wrt b
    Vbf(1:I-1,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))./gridd.dbbayF(:,:,1:2); 
    Vbf(I,:,:)     = (gridd.iiincome(I,:,1:Nz).*(1-param.itax2(I,:,1:Nz)) + r*gridd.bbb(I,:,1:Nz) - param.atax(I,:,1:Nz)).^(-crra);
    % Backward difference
    Vbb(2:I,:,:) = (V(2:I,:,:)-V(1:I-1,:,:))./gridd.dbbayF(:,:,1:2);
    Vbb(1,:,:)   = (gridd.iiincome(1,:,1:Nz).*(1-param.itax2(1,:,1:Nz)) + (guess.r+param.kappa*(gridd.bbb(1,:,1:Nz)<0)).*gridd.bbb(1,:,1:Nz) - param.atax(1,:,1:Nz)).^(-crra); % Boundary condition

    Vbf = max(Vbf,10^(-6));                                                % For numerical stability 
    Vbb = max(Vbb,10^(-6));

    % Consumption and savings with forward difference
    cf = Vbf.^(-1/param.s); 
    sf = gridd.iiincome(:,:,1:2).*(1-param.itax2(:,:,1:2)) + (guess.r+param.kappa*(gridd.bbb(:,:,1:Nz)<0)).*gridd.bbb(:,:,1:Nz) - cf - param.atax(:,:,1:2);
    Hf = (1-param.s)^(-1)*(cf.^(1-param.s)-constant) + Vbf.*sf;
    % Consumption and savings with backward difference
    cb = Vbb.^(-1/param.s); 
    sb = gridd.iiincome(:,:,1:2).*(1-param.itax2(:,:,1:2)) + (guess.r+param.kappa*(gridd.bbb(:,:,1:Nz)<0)).*gridd.bbb(:,:,1:Nz) - cb - param.atax(:,:,1:2);
    Hb = (1-param.s)^(-1)*(cb.^(1-param.s)-constant) + Vbb.*sb;
    % Consumption and derivative of value function at steady state
    c0 = gridd.iiincome(:,:,1:2).*(1-param.itax2(:,:,1:2)) + (guess.r+param.kappa*(gridd.bbb(:,:,1:Nz)<0)).*gridd.bbb(:,:,1:Nz) - param.atax(:,:,1:2);
%     dV0 = c0.^(-param.s);

    if param.AssumeConcavity == 1
        Ib = (sb < 0); 
        If = (sf > 0).*(1-Ib); 
        I0 = (1-If-Ib);
    else 
        I0 = (1-(sf>0)) .* (1-(sb<0)); 
        Iunique = (sb<0).*(1-(sf>0)) + (1-(sb<0)).*(sf>0); 
        Iboth = (sb<0).*(sf>0);
        Ib = Iunique.*(sb<0) + Iboth.*(Hb>Hf);                                    
        If = Iunique.*(sf>0) + Iboth.*(Hb<=Hf);
    end

    c        = cf.*If + cb.*Ib + c0.*I0;
    u        = (1-param.s)^(-1)*(c.^(1-param.s)-param.constant);
    X        = -Ib.*sb./repmat(gridd.dbgridb,1,J,gridd.Nz); % lower
    Y        = -If.*sf./repmat(gridd.dbgridf,1,J,gridd.Nz)+ Ib.*sb./repmat(gridd.dbgridb,1,J,gridd.Nz); % center
    Z        = If.*sf./repmat(gridd.dbgridf,1,J,gridd.Nz);  % upper
    updiag   = zeros(I*J,Nz); 
    lowdiag  = zeros(I*J,Nz);
    centdiag = zeros(I*J,Nz);
    for i = 1:Nz
        centdiag(:,i) = reshape(Y(:,:,i),I*J,1);
    end
    lowdiag(1:I-1,:) = X(2:I,1,:);
    updiag(2:I,:)    = Z(1:I-1,1,:); 
    for j = 2:J
        lowdiag(1:j*I,:) = [lowdiag(1:(j-1)*I,:);squeeze(X(2:I,j,:));zeros(1,Nz)];                                                                    
        updiag(1:j*I,:) = [updiag(1:(j-1)*I,:);zeros(1,Nz);squeeze(Z(1:I-1,j,:))];
    end

    AA1          = spdiags(centdiag(:,1),0,I*J,I*J)+spdiags(updiag(:,1),1,I*J,I*J)+spdiags(lowdiag(:,1),-1,I*J,I*J); 
    AA2          = spdiags(centdiag(:,2),0,I*J,I*J)+spdiags(updiag(:,2),1,I*J,I*J)+spdiags(lowdiag(:,2),-1,I*J,I*J);
    AA           = [AA1, sparse(I*J,I*J); sparse(I*J,I*J), AA2];
    A            = AA + Bswitch; 
    B            = (1/num.Delta + rho + param.kappa)*speye(I*J*Nz) - A;
    u_stacked    = [reshape(u(:,:,1),I*J,1);reshape(u(:,:,2),I*J,1)]; 
    V_stacked    = [reshape(V(:,:,1),I*J,1);reshape(V(:,:,2),I*J,1)];
    vec          = u_stacked + V_stacked/num.Delta;
    V_stacked_12 = B\vec;
    V(:,:,1)     = reshape(V_stacked_12(1:I*J),I,J); 
    V(:,:,2)     = reshape(V_stacked_12(I*J+1:I*J*Nz),I,J);
    Vchange      = V - v;
    v            = V;
    dist(n)      = max(max(max(abs(Vchange))));
    if dist(n)<num.tol_hjb
        break
    end
    end