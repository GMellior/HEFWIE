clear
clc

%% US - undergraduate enrolment

% https://nces.ed.gov/programs/digest/d20/tables/dt20_303.45.asp

opts                       = spreadsheetImportOptions("NumVariables", 9);
opts.Sheet                 = "Digest 2020 Table 303.45";
opts.DataRange             = "A44:I53";
opts.VariableNames         = ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8", "Undergraduate"];
opts.SelectedVariableNames = ["VarName1", "Undergraduate"];
opts.VariableTypes         = ["categorical", "char", "char", "char", "char", "char", "char", "char", "double"];
opts                       = setvaropts(opts, ["Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8"], "WhitespaceRule", "preserve");
opts                       = setvaropts(opts, ["VarName1", "Var2", "Var3", "Var4", "Var5", "Var6", "Var7", "Var8"], "EmptyFieldRule", "auto");
tabn303                    = readtable("tabn303.45.xls", opts, "UseExcel", false);
clear opts

values = tabn303{:, 2}';
v2     = [sum(values(1:2)),values(3:6),sum(values(7:end))];
Age    = {'19 and below','20 and 21','22 to 24','25 to 29','30 to 34','35 years +'};

figure
set(gcf,'Position',[680,144,1140,734])
bar(v2,'FaceColor',[0 0 0.6],'FaceAlpha',0.6,'EdgeColor','none')
title('Enrolment by age 2019')
ylabel('\% of total')
set(gca,'Xtick',1:numel(Age),'XTickLabel',Age)
xtickangle(45)
set(gcf,'Renderer','opengl');
% print(gcf,'EnrollementAgeUnderG2019.eps','-depsc','-r300')

%% UK

% https://www.hesa.ac.uk/news/08-08-2024/sb269-higher-education-student-statistics/numbers

Age    = {'20 and under','21-24','25-29','30 and over'};
years  = {'2018/19','2019/20','2020/21','2021/22','2022/23'};
colmat = flipud([0 0 0.6; 0 0.4 0.8; 0 0.5 1;0.2 0.7 1]); 

data   = [
          71, 70, 67, 68, 68; % 20 and under
          12, 12, 12, 12, 11; % 21-24 years
          6,  6,  7,  7,  6;  % 25-29 years
          11, 12, 14, 14, 15  % 30 years and over
         ];

data = [data(1:2,:);sum(data(3:4,:),1)];
Age2 = {'20 and under','21-24','25 and over'};
figure
set(gcf,'Position',[680,144,1140,734])
h    = bar(data','grouped');
for i = 1:length(h)
    set(h(i),'FaceColor',colmat(i,:),'FaceAlpha',0.6,'EdgeColor','none');
end
title('Enrolment by age (UK) - Entrants and first degree')
ylabel('\% of total')
set(gca,'Xtick',1:numel(years),'XTickLabel',years)
xtickangle(45)
legend(Age2,'Box','Off','Location','northeast')
set(gcf,'Renderer','opengl');
ylim([0 73])
% print(gcf,'EnrollementAgeUKfefd.eps','-depsc','-r300')


%% Average age of enrollment in the OECD (first degree)

opts                  = delimitedTextImportOptions("NumVariables", 2);
opts.DataLines        = [4, Inf];                                          
opts.Delimiter        = ",";
opts.VariableNames    = ["Country", "Age"];
opts.VariableTypes    = ["string", "double"];                              
opts.ExtraColumnsRule = "ignore";
opts.EmptyLineRule    = "read";
opts                  = setvaropts(opts, "Country", "EmptyFieldRule", "auto"); 
OECDAgeofenrolment    = readtable("OECD Age of enrolment.csv", opts);
countries             = cellstr(OECDAgeofenrolment.Country);               
ages                  = OECDAgeofenrolment.Age;                            
countries(4)          = {'Turkey'};
clear opts

figure
set(gcf, 'Position',[680,144,1140,734])
h = bar(ages);
set(h, 'FaceColor', [0 0 0.6], 'FaceAlpha', 0.6, 'EdgeColor', 'none');
title('Average age of enrolment in 2021')
ylabel('Age')
set(gca,'Xtick',1:numel(countries),'XTickLabel',countries)
xtickangle(45)
ylim([17 26])
set(gcf,'Renderer','opengl');
% print(gcf,'EnrollementAgeOECD.eps','-depsc','-r300')

%% Rest

filename                            = 'Total expenditure on educational institutions per full-time equivalent student (2019).xlsx';
sheet                               = 'Sheet1';
data                                = readtable(filename,'Sheet',sheet);
country                             = data.Country; 
country(strcmp(country, 'Türkiye')) = {'Turkey'};
all_tertiary                        = data.AllTertiary; 
all_tertiary_excluding_RnD          = data.ATExcludingR_D; 
growth1219                          = data.Growth1219;
valid_data                          = ~isnan(all_tertiary);
all_tertiary_clean                  = all_tertiary(valid_data);
ateRnD_clean                        = all_tertiary_excluding_RnD(valid_data);
country_clean                       = country(valid_data);
growth                              = growth1219(valid_data);
% Sort data in descending order of all_tertiary
[all_tertiary_sorted,sort_idx]      = sort(all_tertiary_clean,'descend');
country_sorted                      = country_clean(sort_idx);
ATERDsorted                         = ateRnD_clean(sort_idx);
growthsorted                        = growth(sort_idx);

EURcountries = {'Luxembourg','Germany','Austria','Belgium','Netherlands','France', ...
                'Italy','Spain','Portugal','Slovak Republic','Slovenia','Finland', ...
                'Czech Republic','Hungary','Poland','Estonia','Latvia','Lithuania', ...
                'Greece','Ireland','Sweden', 'Norway','Iceland','Denmark','United Kingdom','Turkey'}; 
EUR_indices  = ismember(country_sorted, EURcountries);
EUR_avg      = mean(all_tertiary_sorted(EUR_indices));
EUR_avg_grow = mean(growthsorted(EUR_indices),"omitnan");


figure
set(gcf,'Position',[680,144,1140,734])
bar(all_tertiary_sorted,'FaceColor',[0 0 0.6],'FaceAlpha',0.6,'EdgeColor', 'none')
hold on
bar(ATERDsorted, 'FaceColor', [1 0 0], 'FaceAlpha', 0.6, 'EdgeColor', 'none')
title('Total expenditure per tertiary student in 2019')
ylabel('Expenditure (USD)')
set(gca,'Xtick',1:numel(country_sorted),'XTickLabel',country_sorted)
xtickangle(90)
ax                       = gca;
ax.YAxis.Exponent        = 0;            
ax.YAxis.TickLabelFormat = '%,g';        
yline(EUR_avg, 'k-.', 'LineWidth', 1.5); 
legend('All','Excluding R\&D','All tertirary Europe','Box','Off','Fontsize',14)
hold off
% set(gcf, 'Renderer', 'painters');
% saveas(gcf, 'TertiaryExpStudent2019', 'epsc');
set(gcf,'Renderer','opengl');
% print(gcf,'TertiaryExpStudent2019.eps','-depsc','-r300')

figure;
set(gcf, 'Position',[680,144,1140,734])
bar(growthsorted,'FaceColor',[0 0 0.6],'FaceAlpha',0.6,'EdgeColor', 'none')
hold on
title('Average annual growth 2012-2019 of total expenditure per tertiary student')
ylabel('\% change')
set(gca,'Xtick',1:numel(country_sorted),'XTickLabel',country_sorted)
xtickangle(90)
ax                = gca;
ax.YAxis.Exponent = 0;                   % No comma, no scientific notation
yline(EUR_avg_grow, 'k-.', 'LineWidth', 1.5); 
legend('Country','Europe','Box','Off','Fontsize',14)
hold off
% set(gcf,'Renderer', 'painters');
% saveas(gcf,'TertiaryExpStudentGrowth20122019', 'epsc');
set(gcf,'Renderer','opengl');
% print(gcf,'TertiaryExpStudentGrowth20122019.eps','-depsc','-r300')

%% Scatter plots

filename   = 'OECDfinalCworkingage.csv';
delimiter  = ',';
startRow   = 2;
formatSpec = '%f%C%f%f%f%f%f%[^\n\r]';
fileID     = fopen(filename,'r');
dataArray  = textscan(fileID, formatSpec,'Delimiter',delimiter,'TextType','string','HeaderLines',startRow-1,'ReturnOnError',false,'EndOfLine','\r\n');
fclose(fileID);
datamat    = table(dataArray{1:end-1},'VariableNames',{'VarName1','Country','Year','Public','Private','GINIpost','GINIpre'});
clearvars filename delimiter startRow formatSpec fileID dataArray ans;

latexfont= 1;
public   = table2array(datamat(:,4));
private  = table2array(datamat(:,5));
HEexp    = [public,private];
ginipost = table2array(datamat(:,6));
ginipre  = table2array(datamat(:,7));
fsz      = 20;

%% Post
for k=1:2
    R      = corrcoef(ginipost,HEexp(:,k));
    Rsq(k) = R(1,2);
end

%% Pre
for k=1:2
    R      = corrcoef(ginipre,HEexp(:,k));
    Rsq2(k)= R(1,2);
end

fig1 = figure;
set(fig1,'position',[400 200 1200 500])
subplot(221)
s = scatter(HEexp(:,2),ginipost);
h = lsline;
set(h, 'LineWidth',2,'color',[1 0 0])
s.LineWidth = 1;
s.MarkerEdgeColor = [0 0 0.6];
s.MarkerFaceColor = [0 0 1];
s.MarkerFaceAlpha = 0.5;
ylabel('Gini (post)','Interpreter','Latex','Fontsize',fsz)
if latexfont==1
    set(gca,'Ticklabelinterpreter','Latex','Fontsize',fsz)
else
    set(gca,'Fontsize',fsz)
end
title(sprintf('$\\rho$ = %.2f',Rsq(2)),'Interpreter','Latex','Fontsize',fsz)
grid on
xlim([min(HEexp(:,2)) max(HEexp(:,2))])
ylim([0.23 0.4])

subplot(222)
s = scatter(HEexp(:,1),ginipost);
h = lsline;
set(h, 'LineWidth',2,'color',[1 0 0])
s.LineWidth = 1;
s.MarkerEdgeColor = [0 0 0.6];
s.MarkerFaceColor = [0 0 1];
s.MarkerFaceAlpha = 0.5;
if latexfont==1
    set(gca,'Ticklabelinterpreter','Latex','Fontsize',fsz)
else
    set(gca,'Fontsize',fsz)
end
grid on
title(sprintf('$\\rho$ = %.2f',Rsq(1)),'Interpreter','Latex','Fontsize',fsz)
xlim([min(HEexp(:,1)) max(HEexp(:,1))])
ylim([0.23 0.4])

subplot(223)
s = scatter(HEexp(:,2),ginipre);
h = lsline;
set(h, 'LineWidth',2,'color',[1 0 0])
s.LineWidth = 1;
s.MarkerEdgeColor = [0 0 0.6];
s.MarkerFaceColor = [0 0 1];
s.MarkerFaceAlpha = 0.5;
ylabel('Gini (pre)','Interpreter','Latex','Fontsize',fsz)
xlabel('Private HEexp/GDP','Interpreter','Latex','Fontsize',fsz)
if latexfont==1
    set(gca,'Ticklabelinterpreter','Latex','Fontsize',fsz)
else
    set(gca,'Fontsize',fsz)
end
title(sprintf('$\\rho$ = %.2f',Rsq2(2)),'Interpreter','Latex','Fontsize',fsz)
grid on
xlim([min(HEexp(:,2)) max(HEexp(:,2))])
ylim([0.3 0.57])

subplot(224)
s = scatter(HEexp(:,1),ginipre);
h = lsline;
set(h, 'LineWidth',2,'color',[1 0 0])
s.LineWidth = 1;
s.MarkerEdgeColor = [0 0 0.6];
s.MarkerFaceColor = [0 0 1];
s.MarkerFaceAlpha = 0.5;
xlabel('Public HEexp/GDP','Interpreter','Latex','Fontsize',fsz)
if latexfont==1
    set(gca,'Ticklabelinterpreter','Latex','Fontsize',fsz)
else
    set(gca,'Fontsize',fsz)
end
grid on
title(sprintf('$\\rho$ = %.2f',Rsq2(1)),'Interpreter','Latex','Fontsize',fsz)
xlim([min(HEexp(:,1)) max(HEexp(:,1))])
ylim([0.3 0.57])

%% Another perspective of HE public-private expenditure

fsz  = 20;
fig2 = figure;
set(fig2,'position',[400 200 1200 500])
subplot(121)
s = scatter(HEexp(:,2)./HEexp(:,1),ginipre);
h = lsline;
set(h, 'LineWidth',2,'color',[1 0 0])
s.LineWidth = 1;
s.MarkerEdgeColor = [0 0 0.6];
s.MarkerFaceColor = [0 0 1];
s.MarkerFaceAlpha = 0.5;
ylabel('Gini pre','Interpreter','Latex','Fontsize',fsz)
xlabel('Private HEexp/Public HEexp','Interpreter','Latex','Fontsize',fsz)
if latexfont==1
    set(gca,'Ticklabelinterpreter','Latex','Fontsize',fsz)
else
    set(gca,'Fontsize',fsz)
end
grid on
xlim([min(HEexp(:,2)) max(HEexp(:,2))])
ylim([0.2 0.57])
% Post
subplot(122)
s = scatter(HEexp(:,2)./HEexp(:,1),ginipost);
h = lsline;
set(h, 'LineWidth',2,'color',[1 0 0])
s.LineWidth = 1;
s.MarkerEdgeColor = [0 0 0.6];
s.MarkerFaceColor = [0 0 1];
s.MarkerFaceAlpha = 0.5;
ylabel('Gini post','Interpreter','Latex','Fontsize',fsz)
xlabel('Private HEexp/Public HEexp','Interpreter','Latex','Fontsize',fsz)
if latexfont==1
    set(gca,'Ticklabelinterpreter','Latex','Fontsize',fsz)
else
    set(gca,'Fontsize',fsz)
end
grid on
xlim([min(HEexp(:,2)) max(HEexp(:,2))])
ylim([0.2 0.57])


% Plot 
fsz  = 20;
fig3 = figure;
set(fig3,'position',[400 200 1200 500])
s = scatter(HEexp(:,2)./HEexp(:,1),ginipre);
h = lsline;
set(h, 'LineWidth',2,'color',[1 0 0])
s.LineWidth = 1;
s.MarkerEdgeColor = [0 0 0.6];
s.MarkerFaceColor = [0 0 1];
s.MarkerFaceAlpha = 0.5;
ylabel('Gini','Interpreter','Latex','Fontsize',fsz)
xlabel('Private HEexp/Public HEexp','Interpreter','Latex','Fontsize',fsz)
if latexfont==1
    set(gca,'Ticklabelinterpreter','Latex','Fontsize',fsz)
else
    set(gca,'Fontsize',fsz)
end
grid on
xlim([min(HEexp(:,2)) max(HEexp(:,2))])
ylim([0.3 0.57])

%% Inflation 

load('inflation')
fst = 20;
figure
hold on
plot(d{2}.Dates,d{2}.Values,'linewidth',3,'color',[31/255 19/255 134/255]); % Healthcare
plot(d{1}.Dates,d{1}.Values,'linewidth',3,'color',[0 0 0]);                 % CPI
plot(d{3}.Dates,d{3}.Values,'linewidth',3,'color',[1 0 0]);                 % Tuition
plot(d{4}.Dates,d{4}.Values,'linewidth',3,'color',[18/255 188/255 222/255]);% Housing
datetick('x','yyyy','keeplimits')
recessionplot
hold off
leg = legend('Healthcare','CPI','Tuition','Housing','Location','northwest');
set(leg,'Interpreter','latex','FontSize',24);
legend('boxoff')
ylabel('1982=100','Interpreter','latex')
alpha(0.15)
xlim([min(d{2}.Dates) max(d{2}.Dates)])
set(gca,'FontSize',fst)
set(gca,'XTick', datenum(1970:5:2020,1,1))
set(gca,'XTickLabel', {'1970','1975','1980','1985','1990','1995','2000','2005','2010','2015','2020'});
grid on
