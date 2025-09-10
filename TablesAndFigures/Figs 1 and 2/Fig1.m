clear all
clc
% Create Figure 1
% Set plot properties once
if ~strcmp(get(groot,'defaulttextInterpreter'),'latex')
    set(groot,'defaulttextInterpreter','latex');
    set(groot,'defaultAxesTickLabelInterpreter','latex');
    set(groot,'defaultAxesFontsize',14);
    set(groot,'defaultLegendInterpreter','latex');
    set(groot,'defaultAxesXGrid', 'on');
    set(groot,'defaultAxesYGrid', 'on');
    set(groot,'defaultLineLineWidth',2)
end

%% Aggregate student loans over disposable personal income
load('aggregatesloans');
str3       = 'Federal Student Debt as a Percent of Disposable Personal Income';
%% Distribution of student debt
filename    = 'studentloan.csv';
delimiter   = ',';
formatSpec  = '%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%f%[^\n\r]';
fileID      = fopen(filename,'r');
dataArray   = textscan(fileID, formatSpec,'Delimiter',delimiter,'TextType','string','ReturnOnError',false);
fclose(fileID);
studentloan = table(dataArray{1:end-1}, 'VariableNames', {'VarName1','VarName2','VarName3','VarName4','VarName5','VarName6','VarName7','VarName8','VarName9','VarName10','VarName11','VarName12','VarName13','VarName14','VarName15','VarName16','VarName17','VarName18'});
% Clear temporary variables
clearvars filename delimiter formatSpec fileID dataArray ans;
%% Plots
data        = table2array(studentloan);
Population0 = data(:,[2:2:end]);
Population  = sum(data(:,[2:2:end]),2);
Debt        = sum(data(:,[1:2:end-1]),2);
PopWeights  = data(:,[2:2:end])./repmat(Population,1,size(data,2)/2);
DebtWeights = data(:,[1:2:end-1])./repmat(Debt,1,size(data,2)/2);
MeanDebt    = sum(Debt.*DebtWeights,2)./Population;
fst         = 17; % Ticksize
fsl         = 20; % Fontsize
colors      = [1 0 0; 0 0 0.6];
Quarter     = [2 3 4 1 2 3 4 1 2 3 4 1 2];
Year        = [2017*ones(1,3) 2018*ones(1,4) 2019*ones(1,4) 2020 2020];

X          = categorical({'$<$ 5K','5K to 10K','10K to 20K','20K to 40K','40K to 60K','60K to 80K','80K to 100K','100K to 200K','$>$ 200K'});
X          = reordercats(X,{'$<$ 5K','5K to 10K','10K to 20K','20K to 40K','40K to 60K','60K to 80K','80K to 100K','100K to 200K','$>$ 200K'});
X2         = categorical({'$<$ 10K','10K to 25K','25K to 50K','50K to 100K','100K to 150K','150K to 200K','$>$ 200K'});
X2         = reordercats(X2,{'$<$ 10K','10K to 25K','25K to 50K','50K to 100K','100K to 150K','150K to 200K','$>$ 200K'});
edges      = [0 5 10 20 40 60 80 100 200 250]; 
binIdx      = find(edges(1:end-1)<=MeanDebt(end) & edges(2:end)>MeanDebt(end));
frac        = (MeanDebt(end) - edges(binIdx))/(edges(binIdx+1) - edges(binIdx));
Xnum        = 1:numel(X);
Xpos        = Xnum(binIdx) + frac;

%%
fig = figure('Position',[334,388,1500,420]);
subplot(121)
bar(X,PopWeights(end,:)*100,'edgecolor','none','facecolor',[0 0 0.6],'facealpha',0.6);
hold on
line([Xpos-1/2 Xpos-1/2],ylim,'LineStyle','--','Color','k','LineWidth',1.5);
set(gca,'Fontsize',fst)
set(gca,'XTickLabelRotation',-30)
yl = ylim;
text(binIdx+1,yl(end)*0.95, sprintf('Mean %.2fK',MeanDebt(end)), ...
    'Interpreter', 'latex', ...
    'FontSize', 14, ...
    'VerticalAlignment', 'top', ...
    'HorizontalAlignment', 'left', ...
    'FontWeight', 'normal');
hold off
subplot2=subplot(122);
plot(dates,SDoPI*100,'linewidth',3,'color',[31/255 19/255 134/255])
xlim([728000 737791])
ylim([0 8.1])
datetick('x','yyyy','keeplimits')
recessionplot
set(gca,'Fontsize',fst)
set(subplot2,'FontSize',17,'XTick',...
    [728660 730486 732313 734139 735965 737791],'XTickLabel',...
    {'1995','2000','2005','2010','2015','2020'});
hold off