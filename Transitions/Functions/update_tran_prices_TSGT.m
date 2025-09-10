num.Delta         = dtsvec(t);
guess.t           = t;
guess.r           = r_t(t);
guess.r_t         = r_t;
r                 = r_t(t);
gridd.income      = [wL_t(t)*reprate,wL_t(t),wL_t(t)*zstn,wH_t(t)*reprate,wH_t(t)]; 
switch GT
    case 0
        itax              = itax_t(t);
        param.itax        = [0 itax_t(t) 0 0 itax_t(t)];
    case 1
        itax              = itax_t(t,:);
        param.itax        = [0 itax_t(t,1) 0 0 itax_t(t,2)]; 
end
for jjj=1:NzNed
    gridd.iiincome(:,jjj) = gridd.income(1,jjj)*ones(gridd.I,gridd.J,1);
    param.itax2(:,jjj)    = kron(param.itax(1,jjj),ones(gridd.I,J));
end
P            = P_t(t);
P0           = (1-subsidy)*P;
i_buy        = min(find(gridd.b+abs(gridd.b(1))>P0)); %P0 costs i_buy grid points
i_buymesh2   = gridd.b+abs(gridd.b(1))>P0;
x            = gridd.b-P0;
buildOmega;