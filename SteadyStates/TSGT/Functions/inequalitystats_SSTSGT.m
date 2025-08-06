% Other statistics
Meanb     = (gridd.b.*gridd.bdelta)'*sum(sum(g,3),2);
Varb      = (gridd.b.^2.*gridd.bdelta)'*sum(sum(g,3),2) - Meanb^2;
Varc      = gridd.bdelta'*(sum(sum(g.*(cons.^2),3),2)) - AggC^2;
CVw       = sqrt(Varb)/Meanb;
CVc       = sqrt(Varc)/AggC;
pop       = cumsum(gridd.bdelta.*(sum(sum(g,3),2)));
Lorenzw   = (gridd.b.*gridd.bdelta).*(sum(sum(g,3),2));
Lorenzc   = gridd.bdelta.*(sum(sum(g.*cons,2),3));
Lorenzw   = cumsum(Lorenzw)/sum(Lorenzw);
Lorenzc   = cumsum(Lorenzc)/sum(Lorenzc);
aind      = max(find(gridd.b<0));
% Share of people in debt
debtors   = gridd.bdelta(1:aind,1)'*(sum(sum(g(1:aind,:,:),2),3));
Bg        = trapz(pop,Lorenzw);
Ag        = 0.5-Bg;
Giniw     = Ag/(Ag+Bg);
Bg        = trapz(pop,Lorenzc);
Ag        = 0.5-Bg;
Ginic     = Ag/(Ag+Bg);