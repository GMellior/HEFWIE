% This file builds the Omega matrix for one type
bP                = gridd.b-P0;
affordable        = bP>gridd.bmin;
adjpoints0        = pointslist0(affordable);
rebalancing_list0 = adjpoints0;
sp_addend0        = sparse(I,I);
for mcount = 1:length(rebalancing_list0)
    m = rebalancing_list0(mcount);
    bprime       = max(bP(m),gridd.b(1,1));
    bprime_left  = discretize(bprime, gridd.b);
    bprime_right = bprime_left + 1;
    point11      = pointsmat0(bprime_left,1);
    point12      = pointsmat0(bprime_right,1);
    neighpoints  = [point11 point12];
    totarea    = (gridd.b(bprime_right) - gridd.b(bprime_left));
    weights    = totarea^(-1) * [gridd.b(bprime_right) - bprime bprime - gridd.b(bprime_left)];
    sp_addend0  = sp_addend0 + sparse((m) * ones(2, 1), 0+neighpoints(:), weights(:), I, I);
end
Omega = sp_addend0 + spdiags(~affordable,0,I,I);