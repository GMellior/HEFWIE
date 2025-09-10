
for oloop=1:1
    lext                 = kron([param.lx3;lx1;lx2],ones(I,1));
    vne(:,1)             = [kron(ones(Nz,1),reshape(Vall(:,1,sss),M,1));reshape(Vall(:,2,sss),M,1)];
    vne(1:I,2)           = reshape(bnodebt*((1-parent_abs_small(1,sss))*Vall(:,1,sss)+parent_abs_small(1,sss)*Vall(:,1,mod(sss,2)+1)),M,1);
    vne(1+I:end,2)       = repmat(reshape(bnodebt*((1-parent_abs_small(2,sss))*Vall(:,1,sss)+parent_abs_small(2,sss)*Vall(:,1,mod(sss,2)+1)),M,1),gridd.Nz,1);
    v0pack2              = TSGT_Vedu_transition(param,guess,num,gridd,Aswitch2,vne,lext);
    Vall(:,Nz+1:end,sss) = v0pack2.V;
    V                    = Vall(:,1:Nz,sss);
    VOS_stacked          = V(:);                                           % Fully implicit
    Vold                 = ones(MMu,1);
    for mmmm=1:20
        V = Vall(:,1:Nz,sss);
        % Forward difference for the derivative of v wrt b
        Vbf(1:I-1,:) = (V(2:I,:)-V(1:I-1,:))./gridd.dbyF(:,1:Nz); 
        Vbf(I,:)     = (gridd.iiincome(I,1:Nz).*(1-param.itax2(I,1:Nz)) + guess.r*gridd.bmax).^(-crra); % Boundary condition
        % Backward difference
        Vbb(2:I,:)   = (V(2:I,:)-V(1:I-1,:))./gridd.dbyB(:,1:Nz);
        Vbb(1,:)     = (gridd.iiincome(1,1:Nz).*(1-param.itax2(1,1:Nz)) + (guess.r+param.kappa*(gridd.b(1)<0))*gridd.b(1)).^(-crra); % Boundary condition
    
        Vbf = max(Vbf,10^(-6));                                            % For stability 
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
        indicsubs1 = (c0(:,1)<10^(-6));
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
        X         = -Ib.*sb./repmat(gridd.dbgridb,1,gridd.Nz);             % lower
        Y         = -If.*sf./repmat(gridd.dbgridf,1,gridd.Nz) + Ib.*sb./repmat(gridd.dbgridb,1,gridd.Nz); % center
        Z         = If.*sf./repmat(gridd.dbgridf,1,gridd.Nz);              % upper
        A1        = spdiags(Y(:,1),0,I,I)+spdiags(X(2:I,1),-1,I,I)+spdiags([0;Z(1:I-1,1)],1,I,I);
        A2        = spdiags(Y(:,2),0,I,I)+spdiags(X(2:I,2),-1,I,I)+spdiags([0;Z(1:I-1,2)],1,I,I);
        A         = [A1,sparse(I,I);sparse(I,I),A2] + Aswitch1;
        B         = (1/num.Delta + rho + param.kappa)*speye(I*J*Nz) - A;
        u_stacked = [reshape(u(:,1),I,1);reshape(u(:,2),I,1)];
        V_stacked = [reshape(V(:,1),I,1);reshape(V(:,2),I,1)];
    
        Vsaux1              = Omega*v0pack2.V(:,1);
        Vsaux1(1:i_buy-1,1) = Vsaux1(1:i_buy-1,1)-5000;
        Vsaux2              = Vsaux1;
        Vstar_stacked       = [Vsaux1;Vsaux2] - param.mcost(1);
        vec                 = u_stacked + VOS_stacked/num.Delta + param.altruism*param.kappa*repmat(bnodebt*( (1-parent_abs_small(1,sss))*Vall(:,1,sss)+ parent_abs_small(1,sss)*Vall(:,1,mod(sss,2)+1)),gridd.Nz,1)+0.001*max(Vstar_stacked-VOS_stacked,0);
        q                   = -vec + B*Vstar_stacked; 
        z0                  = V_stacked-Vstar_stacked;
        z                   = LCP(B,q,lbnd,ubnd,z0,0);
        LCP_error           = max(abs(z.*(B*z + q)));
        V_stacked           = z+Vstar_stacked; 
        V                   = reshape(V_stacked,I,Nz);
        vdiff               = max(abs(V_stacked-Vold));
        Vold                = V_stacked;
        
        if vdiff<num.tol_hjb
            break
        end
    end
    Vall(:,1:Nz,sss) = V;
end
loopcount(:,:,t) = [mmmm,vdiff;v0pack2.nloop,v0pack2.dist;oloop,NaN];