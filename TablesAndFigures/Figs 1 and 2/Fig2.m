
clear all
clc

loadOECD;
countries   = OECD.country';
delrow = [];

%% Drop unused countries
ctrydel = ["BRA","CRI","ISR","IND","IDN","RUS","SAU","ZAF","LUX","ARG","ISL","HUN","LTU","LTA","EST","MEX",...
            "TUR","LVA","SVK","SVN","GRC","COL","CHL","CZE","POL","IRL","PRT","BEL","AUS","NZL","CHN","CHE"];
for j=1:numel(ctrydel)
    for i=1:numel(OECD.country)
        if ctrydel(j)==countries(i)
                delrow = [delrow i];
        end
    end
end

OECD(delrow,:) = [];

countries   = OECD.country';
ctryunique  = unique(countries,'stable');
yearunique  = unique(OECD.Year);
data        = nan(numel(yearunique),1+numel(ctryunique));
data(:,1)   = yearunique;

for j=1:numel(ctryunique)
    for i=1:numel(OECD.country)
        if ctryunique(j)==countries(i)
            for k=1:4
                if OECD.Year(i)==yearunique(k)
                    data(k,j+1) = OECD.Value(i);
                end
            end
        end
    end
end

ds2013 = [ctryunique;data(1,2:end)]';
ds2013 = sortrows(ds2013,2);
out2013 = ds2013(all(~isnan(double(ds2013(:,2))),2),:);
ds2014 = [ctryunique;data(2,2:end)]';
ds2014 = sortrows(ds2014,2);
out2014 = ds2014(all(~isnan(double(ds2014(:,2))),2),:);
ds2015 = [ctryunique;data(3,2:end)]';
ds2015 = sortrows(ds2015,2);
out2015 = ds2015(all(~isnan(double(ds2015(:,2))),2),:);

loadOECD2;
countries   = OECDpubpriv.country';
delrow = [];


for j=1:numel(ctrydel)
    for i=1:numel(OECDpubpriv.country)
        if ctrydel(j)==countries(i)
                delrow = [delrow i];
        end
    end
end

OECDpubpriv(delrow,:) = [];
OECDpubpriv(find(OECDpubpriv.pubpriv=="PRIV"),:) = [];
OECDpubpriv(:,2) = [];

countries   = OECDpubpriv.country';
ctryunique  = unique(countries,'stable');
yearunique  = unique(OECDpubpriv.year);
data        = nan(numel(yearunique),1+numel(ctryunique));
data(:,1)   = yearunique;

for j=1:numel(ctryunique)
    for i=1:numel(OECDpubpriv.country)
        if ctryunique(j)==countries(i)
            for k=1:4
                if OECDpubpriv.year(i)==yearunique(k)
                    data(k,j+1) = OECDpubpriv.value(i);
                end
            end
        end
    end
end

yrdel = [2013,2014,2016]'; % We plot 2015
clear('delrow')
delrow = [];
for j=1:numel(yrdel)
    for i=1:numel(OECDpubpriv.year)
        if yrdel(j)==OECDpubpriv.year(i)
                delrow = [delrow i];
        end
    end
end

OECDpubpriv(delrow,:) = [];
OECD(find(OECD.Year==2013),:) = [];
OECD(find(OECD.Year==2014),:) = [];
OECD(find(OECD.Year==2016),:) = [];
store0 = nan(numel(OECD(:,1)),1);
store1 = nan(numel(OECD(:,1)),1);


for j=1:numel(OECD(:,1))
    for i=1:numel(OECDpubpriv(:,1))
        if OECD.country(j)==OECDpubpriv.country(i)
                store0(j) = OECD.Value(j);
                store1(j) = OECD.Value(j)*OECDpubpriv.value(i)/100;
                store2(j) = OECD.country(j);
        end
    end
end

store0(find(store2=="DNK")) = []; % DNK is missing a value
store1(find(store2=="DNK")) = [];
store2(find(store2=="DNK")) = [];
dataplot = sortrows([store2' store0 store1],2);

fszx      = 18; % Fontsize x axis
fszl      = 18; % Fontsize legend
falph     = 0.7;
latexfont = 1;
fig2 = figure(2);
set(fig2,'position',[400 200 1200 500])
bar(double(dataplot(:,2)),'FaceAlpha',falph,'EdgeColor','none')
% axis tight
ylim([0 2.75])
set(gca,'XTick',1:size(store1,1)) 
set(gca,'XTickLabel',dataplot(:,1));
set(gca,'Fontsize',fszx);
grid on
hold on
bar(double(dataplot(:,3)),'r','FaceAlpha',falph,'EdgeColor','none')
leg=legend('Private','Public');
set(leg,'Interpreter','Latex','Fontsize',fszl,'location','Northwest');
legend('boxoff')
% axis tight
ylim([0 2.75])
set(gca,'XTick',1:size(store1,1)) 
set(gca,'XTickLabel',dataplot(:,1));
set(gca,'Fontsize',fszx);
if latexfont==1
    set(gca,'Ticklabelinterpreter','Latex','Fontsize',fszx)
else
    set(gca,'Fontsize',fszx)
end
grid on
hold off