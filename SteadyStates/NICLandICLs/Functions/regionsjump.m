% Where you jump to after discrete choices

aAdj        = zeros(I,J);
bAdj        = aAdj;
sumGridP    = round(gridd.bb-P0 - max(gridd.aminfed-max(gridd.aa,gridd.aminfed),-P0),7);
sumGridP(sumGridP<gridd.b(1)) = -5000;
bc          = gridd.b-P0 + max(gridd.aa-gridd.aminfed,0);
aAdj        = (gridd.aa-P0).*(gridd.aa-P0>=gridd.aminfed)...
              + (gridd.aa-P0<gridd.aminfed).*( (bc>=gridd.b(1)).*(gridd.aminfed.*(gridd.aa>=gridd.aminfed) + gridd.aa.*(gridd.aa<gridd.aminfed)) + gridd.aa.*(bc<gridd.b(1)));
bAdj        = gridd.bb.*(gridd.aa-P0>=gridd.aminfed)...
                      + (gridd.aa-P0 <gridd.aminfed).*(bc.*(bc>=gridd.b(1)) + gridd.bb.*(bc<gridd.b(1)));
bAdj        = repmat(bAdj,1,1,2);
aAdj        = repmat(aAdj,1,1,2);   