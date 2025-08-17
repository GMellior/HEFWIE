
Meanb     = gridd.da*(gridd.b.*gridd.bdelta)'*(sum(sum(sum(g,4),3),2));
Varb      = gridd.da*((gridd.b.^2).*gridd.bdelta)'*(sum(sum(sum(g,4),3),2)) - Meanb^2;
Vara      = gridd.da*((gridd.aa(:,1).^2).*gridd.bdelta)'*(sum(sum(sum(g,4),3),2)) - studentdebt^2;
Varc      = gridd.da*gridd.bdelta'*(sum(sum(sum(g.*(cons.^2),2),3),4)) - AggC^2;
CVb       = sqrt(Varb)/Meanb;
CVa       = sqrt(Vara)/abs(studentdebt);
CVc       = sqrt(Varc)/AggC;
pop       = cumsum(gridd.da*(gridd.bdelta).*(sum(sum(sum(g,4),3),2)));
Lorenzw   = gridd.da*(gridd.b.*gridd.bdelta).*(sum(sum(sum(g,4),3),2));
Lorenzw2  = gridd.da*(gridd.bdelta).*(sum(sum(sum(g.*repmat(gridd.bb(:,:,1)+gridd.aa(:,:,1),1,1,NzNed,abtypes),3),2),4));
Lorenzc   = gridd.da*(gridd.bdelta).*(sum(sum(sum(g.*cons, 2),3),4));
Lorenzw   = cumsum(Lorenzw)/sum(Lorenzw);
Lorenzw2  = cumsum(Lorenzw2)/sum(Lorenzw2);
Lorenzc   = cumsum(Lorenzc)/sum(Lorenzc);
% Stock of people in debt
aind      = max(find(gridd.b<0));
Bdebtors  = gridd.da*(gridd.bdelta(1:aind,1))'*(sum(sum(sum(g(1:aind,:,:,:),2),3),4));
Adebtors  = gridd.da*(gridd.bdelta)'*(sum(sum(sum(g.*(repmat(gridd.aa,1,1,NzNed,abtypes)<0),3),2),4));
ABdebtors = gridd.da*(gridd.bdelta)'*(sum(sum(sum(g.*(repmat(gridd.aa+gridd.bb,1,1,NzNed,abtypes)<0),3),2),4));


Bg      = trapz(pop,Lorenzw);
Ag      = 0.5-Bg;
Giniw   = Ag/(Ag+Bg);
Bg      = trapz(pop,Lorenzw2);
Ag      = 0.5-Bg;
Giniw2  = Ag/(Ag+Bg);
Bg      = trapz(pop,Lorenzc);
Ag      = 0.5-Bg;
Ginic   = Ag/(Ag+Bg);
