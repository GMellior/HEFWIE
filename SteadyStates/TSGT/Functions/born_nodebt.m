% This file builds the Omega matrix for being born without debt

bpos              = gridd.b>0;
bP                = gridd.b.*bpos + 0*(~bpos);
pointslist0       = (1:I)';
affordable        = ~bpos;
adjpoints0        = pointslist0(affordable);
rebalancing_list0 = adjpoints0;
sp_addend0        = sparse(I,I);
pointsmat0        = reshape(pointslist0, I, 1);
for mcount = 1:length(rebalancing_list0)
    m = rebalancing_list0(mcount);
    bprime        = max(bP(m),0);
    bprime_left   = discretize(bprime, gridd.b); % This gives the index of the left edge for b'
    bprime_right  = bprime_left + 1;
    % Map from grid indexes to points
    point11       = pointsmat0(bprime_left,1);
    point12       = pointsmat0(bprime_right,1);
    neighpoints   = [point11 point12];
    totarea       = (gridd.b(bprime_right) - gridd.b(bprime_left));
    weights       = totarea^(-1) * [gridd.b(bprime_right) - bprime bprime - gridd.b(bprime_left)];
    sp_addend0    = sp_addend0 + sparse((m) * ones(2, 1), 0+neighpoints(:), weights(:), I, I);
end
bnodebt          = sp_addend0 + spdiags(~affordable,0,I,I);