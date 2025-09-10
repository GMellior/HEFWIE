
pop            = cumsum(gridd.bdelta.*(sum(sum(full(g),3),2)));
Lorenzw        = (gridd.b.*gridd.bdelta).*(sum(sum(full(g),3),2));
Lorenzc        = gridd.bdelta.*(sum(sum(full(g).*cons_t(:,:,:,n),2),3));
Lorenzw        = cumsum(Lorenzw)/sum(Lorenzw);
Lorenzc        = cumsum(Lorenzc)/sum(Lorenzc);
debtors_t(n,1) = gridd.bdelta(1:aind,1)'*(sum(sum(full(g(1:aind,:,:)),2),3));
Bg             = trapz(pop,Lorenzw);
Ag             = 0.5-Bg;
Giniw_t(n,1)   = Ag/(Ag+Bg);
Bg             = trapz(pop,Lorenzc);
Ag             = 0.5-Bg;
Ginic_t(n,1)   = Ag/(Ag+Bg);