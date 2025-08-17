
switch Picksystem
    case 0
        aDrift        = param.atax;
        aDrift(:,:,3) = ra*gridd.aa.*(gridd.aa>gridd.aa(1,1))+param.atax(:,:,3); 
        aDrift(:,:,5) = ra*gridd.aa.*(gridd.aa>gridd.aa(1,1))+param.atax(:,:,5);
    case 1
        aDrift        = ra*repmat(gridd.aa,1,1,NzNed)+param.atax;
        aDrift(:,:,3) = unsubsidized*aDrift(:,:,3); 
    case 2
        aDrift  = -param.amort*repmat(gridd.aa,1,1,NzNed);
end
%%%
if max(max(max(ra*gridd.a+param.atax(1:1,1,5),0)))<0
    warning('Educated and employed with amin cannot repay student loans')
end
% The government has to prevent debt going deeper than amin - so they cover interest at \bar{\bar{a}}
% The loop below puts all the aDrifts in the Bswitch matrix - capture student loan flows
BauxUN = [];
for jj=1:NzNed
    chi    = -min(aDrift(:,:,jj),0)/da; %
    yy     = - max(aDrift(:,:,jj),0)/da + min(aDrift(:,:,jj),0)/da;
    zeta   = max(aDrift(:,:,jj),0)/da;
    updiag = zeros(I,1); 
    for j=1:J
        updiag=[updiag;zeta(:,j)];
    end
    centdiag=chi(:,1)+yy(:,1);
    for j=2:J-1
        centdiag=[centdiag;yy(:,j)];
    end
    centdiag = [centdiag;yy(:,J)+zeta(:,J)];
    lowdiag  = chi(:,2);
    for j=3:J
        lowdiag=[lowdiag;chi(:,j)]; 
    end
    if (jj==1)||(jj==4)
        Baux1   = spdiags(centdiag,0,M,M)+spdiags(lowdiag,-I,M,M)+spdiags(updiag,I,M,M);
    elseif (jj==2)||(jj==5)
        Baux2   = spdiags(centdiag,0,M,M)+spdiags(lowdiag,-I,M,M)+spdiags(updiag,I,M,M);
    elseif jj==3
        Bauxstn = spdiags(centdiag,0,M,M)+spdiags(lowdiag,-I,M,M)+spdiags(updiag,I,M,M);
    end
    if jj==2
        Bswitch  = [Baux1 - speye(M)*la(1), speye(M)*la(1);speye(M)*la(2),Baux2 - speye(M)*la(2)];                                       
    elseif jj==5
        Bswitch2 = [Bauxstn - speye(M)*(param.la5+param.la6),speye(M)*param.la5,speye(M)*param.la6;...
                    sparse(M,M), Baux1 - speye(M)*param.la3, speye(M)*param.la3;...
                    sparse(M,M),speye(M)*param.la4,Baux2 - speye(M)*param.la4];
    end
end
Bswitch  = Bswitch  + ICL2switchUn;
Bswitch2 = Bswitch2 + ICL2switchEd;