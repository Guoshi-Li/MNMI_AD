
% Plot network-averaged recurrent excitation and inhibition together for
% NC, MCI and AD
% Plot spontaneous input for the 3 groups

clc;
close all;
clear all;

load PARA_NC_MCI;

W_NC2 = W_NM;
W_MCI = W_MD;

load PARA_MCI_AD;

W_NC = W_NC2;
W_AD = W_MD;

load AD_SC10;
load ADLabelNet2;

n1 = 48;  % number of NC subjects
n2 = 48;  % number of MCI subjects
n3 = 48;  % number of AD subjects

NR = 46;  % number of regions
CN = NP;  % number of inter-regional EC

NP = 2*NR+1+CN;  % Number of total parameters
NM = 4;          % Number of networks
Nb = N4;


%==========================================================================
%              Individual Connection Analysis
%==========================================================================

% Average weight
Wm_NC  = mean(W_NC);
Wm_MCI = mean(W_MCI);
Wm_AD  = mean(W_AD);

STD_NC  = std(W_NC);
STD_MCI = std(W_MCI);
STD_AD  = std(W_AD);

SE_NC  = STD_NC/sqrt(n1);
SE_MCI = STD_MCI/sqrt(n1);
SE_AD  = STD_AD/sqrt(n3);

% Put NC and MDD together
Wm_NC_MD  = [Wm_NC;  Wm_MCI; Wm_AD];
STD_NC_MD = [STD_NC; STD_MCI; STD_AD];
SE_NC_MD  = [SE_NC; SE_MCI; SE_AD];


% Divide the parameters into four sets
Wm_EE = Wm_NC_MD(:, 1:NR);
Wm_IE = Wm_NC_MD(:, NR+1:2*NR);
Wm_GC = Wm_NC_MD(:, 2*NR+1:NP-1);
Wm_PI = Wm_NC_MD(:, NP);


% Standard Error
EE_ER = SE_NC_MD(:, 1:NR);
IE_ER = SE_NC_MD(:, NR+1:2*NR);
GC_ER = SE_NC_MD(:, 2*NR+1:NP-1);


% Re-arrange the ROIs according to the fucntional network

Wm_EE_SORT = [];
Wm_IE_SORT = [];

EE_ER_SORT = [];
IE_ER_SORT = [];

for L = 1:NM
    
    idx = NID{L};
    Wm_EE_SORT = [Wm_EE_SORT Wm_EE(:, idx)];
    Wm_IE_SORT = [Wm_IE_SORT Wm_IE(:, idx)];
    
    EE_ER_SORT = [EE_ER_SORT EE_ER(:, idx)];
    IE_ER_SORT = [IE_ER_SORT IE_ER(:, idx)];

end
       

Wm_EE_SORT = Wm_EE_SORT';
Wm_IE_SORT = Wm_IE_SORT';

EE_ER_SORT = EE_ER_SORT';
IE_ER_SORT = IE_ER_SORT';



% %========================================================================
% %                        Network Analysis
% %========================================================================
% 
WNee = zeros(3, NM);
WNie = zeros(3, NM);

% Intra-regional weights
for i = 1:NM
    
    WTe = [];
    WTi = [];
    
    idx = NID{i};
    WTe = Wm_EE(:,idx);
    WTi = Wm_IE(:,idx);
    
    WNee(:,i) = mean(WTe, 2);
    WNie(:,i) = mean(WTi, 2);
   
end


Wee_NC = W_NC(:,1:NR);
Wie_NC = W_NC(:,NR+1:2*NR);
Wgc_NC = W_NC(:,2*NR+1:NP-1);
SPI_NC = W_NC(:,NP);

Wee_MCI = W_MCI(:,1:NR);
Wie_MCI = W_MCI(:,NR+1:2*NR);
Wgc_MCI = W_MCI(:,2*NR+1:NP-1);
SPI_MCI = W_MCI(:,NP);

Wee_AD = W_AD(:,1:NR);
Wie_AD = W_AD(:,NR+1:2*NR);
Wgc_AD = W_AD(:,2*NR+1:NP-1);
SPI_AD = W_AD(:,NP);

WeeN_NC = zeros(n1, NM);
WieN_NC = zeros(n1, NM);

WeeN_MCI = zeros(n2, NM);
WieN_MCI = zeros(n2, NM);

WeeN_AD = zeros(n3, NM);
WieN_AD = zeros(n3, NM);


for i = 1:NM
    
    idx = NID{i};
    
    % For NC
    
    WTe = [];
    WTi = [];  
       
    WTe = Wee_NC(:,idx);
    WTi = Wie_NC(:,idx);
    
    WeeN_NC(:,i) = mean(WTe, 2);
    WieN_NC(:,i) = mean(WTi, 2);
    
   % For MCI 
    
    WTe = [];
    WTi = [];  
       
    WTe = Wee_MCI(:,idx);
    WTi = Wie_MCI(:,idx);
    
    WeeN_MCI(:,i) = mean(WTe, 2);
    WieN_MCI(:,i) = mean(WTi, 2);


  % For AD 
    
    WTe = [];
    WTi = [];  
       
    WTe = Wee_AD(:,idx);
    WTi = Wie_AD(:,idx);
    
    WeeN_AD(:,i) = mean(WTe, 2);
    WieN_AD(:,i) = mean(WTi, 2);
     
        
end


WN_NC  = [WeeN_NC  WieN_NC  SPI_NC];
WN_MCI = [WeeN_MCI WieN_MCI SPI_MCI];
WN_AD  = [WeeN_AD  WieN_AD  SPI_AD];

WeeN_NC_Mean  = mean(WeeN_NC);
WeeN_MCI_Mean = mean(WeeN_MCI);
WeeN_AD_Mean  = mean(WeeN_AD);

WieN_NC_Mean  = mean(WieN_NC);
WieN_MCI_Mean = mean(WieN_MCI);
WieN_AD_Mean  = mean(WieN_AD);

WeeN = [ WeeN_NC_Mean; WeeN_MCI_Mean; WeeN_AD_Mean];
WieN = [ WieN_NC_Mean; WieN_MCI_Mean; WieN_AD_Mean];


% Standard errors
STD_NET_NC  = std(WN_NC);
STD_NET_MCI = std(WN_MCI);
STD_NET_AD  = std(WN_AD);

SE_NET_NC  = STD_NET_NC/sqrt(n1);
SE_NET_MCI = STD_NET_MCI/sqrt(n2);
SE_NET_AD  = STD_NET_AD/sqrt(n3);

SE_NET_NC_MD = [SE_NET_NC; SE_NET_MCI; SE_NET_AD];

SE_EE_NET = SE_NET_NC_MD(:,1:NM);
SE_IE_NET = SE_NET_NC_MD(:,NM+1:2*NM);
SE_PI = SE_NET_NC_MD(:,end);

%% ========================================================================
% Statistical Test

q = 0.05;

SPI_ALL = [SPI_NC SPI_MCI SPI_AD];

TV = [];
PV = [];

for i = 1:2
   
   V1 = SPI_ALL(:,i);

   for j = (i+1):3

      V2 = SPI_ALL(:,j);
      [h,p,ci,stats1]  = ttest2(V1, V2);
      TV = [TV; stats1.tstat];
      PV = [PV; p];

   end

end
  
TV
PV

[pID1,pN1] = gretna_FDR(PV,q)

FDR1 = mafdr(PV,'BHFDR', 1)


%==========================================================================
%                         Ploting
%==========================================================================

Wm_EE  = Wm_EE';
Wm_IE  = Wm_IE';

Wm_PI = Wm_PI';

WNee = WNee';
WNie = WNie';

WeeN = WeeN';
WieN = WieN';


% STD
EE_ER = EE_ER';
IE_ER = IE_ER';

SE_EE_NET = SE_EE_NET';
SE_IE_NET = SE_IE_NET';
SE_PI = SE_PI';



%% ========================================================================
% Generate bar graphs

ROI_LABEL_NET = {'DMN', 'SAL', 'EXE', 'LIM' };                      

Xn = 1:NM;

figure;
bar(Xn, WeeN);
hold on;
ngroups = size(WeeN, 1);
nbars = size(WeeN, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, WeeN(:,i), SE_EE_NET(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 14);
title('Recurrent Excitation', 'FontSize', 16);
set(gca, 'XTick', Xn, 'XTickLabel', ROI_LABEL_NET);
ylabel('Weight', 'FontSize', 14);
axis([0.5 4.5 2.0 3.5]);
set(gca, 'YTick', [2:0.5:3.5]); 
box('off');
legend('NC', 'MCI', 'AD');
xtickangle(45);


figure;
bar(Xn, WieN);
hold on;
ngroups = size(WieN, 1);
nbars = size(WieN, 2);
% Calculating the width for each bar group
groupwidth = min(0.8, nbars/(nbars + 1.5));
for i = 1:nbars
    x = (1:ngroups) - groupwidth/2 + (2*i-1) * groupwidth / (2*nbars);
    er = errorbar(x, WieN(:,i), SE_IE_NET(:,i), '.');
    er.Color = [0 0 0];
    er.LineWidth = 1;
    er.LineStyle = 'none'; 
end
hold off
set(gca, 'FontSize', 14);
title('Recurrent Inhibition', 'FontSize', 16);
set(gca, 'XTick', Xn, 'XTickLabel', ROI_LABEL_NET);
ylabel('Weight', 'FontSize', 14);
axis([0.5 4.5 2.0 3.5]);
set(gca, 'YTick', [2:0.5:3.5]); 
box('off');
legend('NC', 'MCI', 'AD');
xtickangle(45);



%==========================================================================
% Generate boxplots

DATA_WEE = {WeeN_NC, WeeN_MCI, WeeN_AD};
DATA_WIE = {WieN_NC, WieN_MCI, WieN_AD};

figure;
% boxplotGroup(DATA_WEE)
% h = boxplotGroup(DATA_WEE,'Colors','gbk','GroupType','betweenGroups') % default
h = boxplotGroup(DATA_WEE); % default
axis([0, 16, 2, 4]);
set(gca,'FontSize', 14);
groupCenters = @(nGroups,nMembers,interGroupSpace) ...
    nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;
x1CenterTicks = groupCenters(numel(DATA_WEE), size(DATA_WEE{1},2), 1);
set(gca,'XTick',x1CenterTicks,'XTickLabels',{'DMN','SAL','EXE','LIM'}, 'FontSize', 16)
legend(h.boxchartGroup, 'NC', 'MCI', 'AD', 'FontSize', 16);
ylabel('Weight', 'FontSize', 16);


figure;
% boxplotGroup(DATA_WEE)
% h = boxplotGroup(DATA_WEE,'Colors','gbk','GroupType','betweenGroups') % default
h = boxplotGroup(DATA_WIE); % default
axis([0, 16, 2, 4]);
set(gca,'FontSize', 14);
groupCenters = @(nGroups,nMembers,interGroupSpace) ...
    nGroups/2+.5 : nGroups+interGroupSpace : (nGroups+interGroupSpace)*nMembers-1;
x1CenterTicks = groupCenters(numel(DATA_WIE), size(DATA_WIE{1},2), 1);
set(gca,'XTick',x1CenterTicks,'XTickLabels',{'DMN','SAL','EXE','LIM'}, 'FontSize', 16)
legend(h.boxchartGroup, 'NC', 'MCI', 'AD', 'FontSize', 16);
ylabel('Weight', 'FontSize', 16);


% Spontaneous input
Xn1 = 1:3;

figure;
bar(Xn1, Wm_PI);
hold on;
er = errorbar(Xn1, Wm_PI, SE_PI, '.');
er.Color = [0 0 0];
er.LineWidth = 1;
er.LineStyle = 'none'; 
hold off
set(gca, 'FontSize', 14);
title('Spont. Input', 'FontSize', 16);
set(gca, 'XTick', [1, 2, 3], 'XTickLabel', {'NC', 'MCI', 'AD'});
axis([0.0 4.0 0.0 0.4]);
box('off');
% legend('NC', 'MDD');

