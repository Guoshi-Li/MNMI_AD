
% Plot intra-/inter-network EC

clc;
close all;
clear all;

Flag = 1; % 1: NC vs. MCI; 2: NC vs AD; 3: MCI vs AD

if (Flag == 1)
   load PARA_NC_MCI;
elseif (Flag == 2)
   load PARA_NC_AD;
else
   load PARA_MCI_AD;
end

load AD_SC10;
load ADLabelNet2;

NR = 46;    % Number of regions
CN = NP;    % Number of inter-regional EC

NP = 2*NR+1+CN;  % Number of total parameters
N  = 4;          % Number of networks
Nb = N4;

q = 0.05;


n1 = 48;  % Number of subjects
n2 = 48;  % Number of subjects

GC_NC = zeros(NR, NR, n1);
GC_MD = zeros(NR, NR, n2);

Combined_GC_NC = zeros(N, N, n1);
Combined_GC_MD = zeros(N, N, n2);

EXC_GC_NC = zeros(N, N, n1);
INH_GC_NC = zeros(N, N, n1);

EXC_GC_MD = zeros(N, N, n2);
INH_GC_MD = zeros(N, N, n2);

VE_NC = [];
VI_NC = [];
VC_NC = [];

VE_MD = [];
VI_MD = [];
VC_MD = [];


% For NC/MCI

for k = 1:n1
    m = 1;
    
    TEMP1_GC_NC = [];
    TEMP2_GC_NC = [];
    
    GC_NC_UNSORTED = zeros(NR, NR);
    
    Wgc = W_NM(k, 2*NR+1:NP-1);
    
    for j=1:NR
      for i=1:NR
       
         if (MAP(i,j)==1)
           
            GC_NC_UNSORTED(i,j) = Wgc(m);
            m = m+1;
         end      
        
      end  
    end
    
    
   % Sort the matrix so that ROIs belonging to the same network are next to each other 
    for L = 1:N
    
       idx = NID{L};
       TEMP1_GC_NC = [TEMP1_GC_NC; GC_NC_UNSORTED(idx,:)];
         
    end
    
   
   for L = 1:N
    
      idx = NID{L};
      TEMP2_GC_NC = [TEMP2_GC_NC TEMP1_GC_NC(:, idx)];
          
   end
 
  
% =========================================================================
% Exchange cACC (15-16) with PCC (17-18)

GC0_NC = TEMP2_GC_NC; 

GC1_NC = [GC0_NC(1:14, :); GC0_NC(17:18, :); GC0_NC(15:16, :); GC0_NC(19:end,:)];
GC1_NC = [GC1_NC(:,1:14)   GC1_NC(:,17:18)   GC1_NC(:,15:16)   GC1_NC(:, 19:end)];


% =========================================================================

   GC_NC(:,:,k) = GC1_NC;
   GC = GC1_NC; 
      
    
  % Combine ROIs into different networks
  
  mi1 = 1;
  mj1 = 1;
  
  for i = 1:N
       
     mj1 = 1;
     mi2 = mi1+Nb(i)-1; 
     
    for j = 1:N
     
     SBLOCK =0;    
          
     mj2 = mj1+Nb(j)-1;
              
     BLOCK = GC(mi1:mi2, mj1:mj2);
     
     [a, b]   = find(BLOCK~=0);
     [l1, l2] = find(BLOCK > 0);
     [r1, r2] = find(BLOCK < 0); 
     
     % Combined
     if ~isempty(a)
       Combined_GC_NC(i, j, k) = sum(sum(BLOCK))/length(a);
     end    
     
     % Excitatory 
     Ln1 = length(l1);
     
     if ~isempty(l1)
       for w = 1:Ln1
          SBLOCK = SBLOCK + BLOCK( l1(w), l2(w) );
       end
     
       EXC_GC_NC (i, j, k) = SBLOCK/Ln1;
     end
     
     
     % Inhibitory 
     SBLOCK = 0;
     
     Ln1 = length(r1);
     
     if  ~isempty(r1)
       for w = 1:Ln1
          SBLOCK = SBLOCK + BLOCK( r1(w), r2(w) );
       end
     
       INH_GC_NC (i, j, k) = SBLOCK/Ln1;
     end
     
     
     mj1 = mj1+Nb(j);
    
     
    end   
     
    mi1 = mi1+Nb(i);
    
  end
    
  TEMP0 = Combined_GC_NC(:,:,k); 
  TEMP1 = EXC_GC_NC(:,:,k);
  TEMP2 = INH_GC_NC(:,:,k);   
  
  VC_NC = [VC_NC TEMP0(:)];
  VE_NC = [VE_NC TEMP1(:)];
  VI_NC = [VI_NC TEMP2(:)];
    
    
    
end




% For MCI/AD

for k = 1:n2
    m = 1;
    
    TEMP1_GC_MD = [];
    TEMP2_GC_MD = [];
    
    GC_MD_UNSORTED = zeros(NR, NR);
    
    Wgc = W_MD(k, 2*NR+1:NP-1);
    
    for j=1:NR
      for i=1:NR
       
         if (MAP(i,j)==1)
           
            GC_MD_UNSORTED(i,j) = Wgc(m);
            m = m+1;
         end      
        
      end  
    end
    
    
   % Sort the matrix so that ROIs belonging to the same network are next to each other 
    for L = 1:N
    
       idx = NID{L};
       TEMP1_GC_MD = [TEMP1_GC_MD; GC_MD_UNSORTED(idx,:)];
         
    end
    
   
   for L = 1:N
    
      idx = NID{L};
      TEMP2_GC_MD = [TEMP2_GC_MD TEMP1_GC_MD(:, idx)];
          
   end
 

% =========================================================================
% Exchange cACC (15-16) with PCC (17-18)

GC0_MD = TEMP2_GC_MD; 

GC1_MD = [GC0_MD(1:14, :); GC0_MD(17:18, :); GC0_MD(15:16, :); GC0_MD(19:end,:)];
GC1_MD = [GC1_MD(:,1:14)   GC1_MD(:,17:18)   GC1_MD(:,15:16)   GC1_MD(:, 19:end)];


% =========================================================================

   GC_MD(:,:,k) = GC1_MD;
   GC = GC1_MD; 
      
    
  % Combine ROIs into different networks
  
  mi1 = 1;
  mj1 = 1;
  
  for i = 1:N
       
     mj1 = 1;
     mi2 = mi1+Nb(i)-1; 
     
    for j = 1:N
     
     SBLOCK =0;    
          
     mj2 = mj1+Nb(j)-1;
              
     BLOCK = GC(mi1:mi2, mj1:mj2);
     
     [a, b]   = find(BLOCK~=0);
     [l1, l2] = find(BLOCK > 0);
     [r1, r2] = find(BLOCK < 0); 
     
     % Combined
     if ~isempty(a)
       Combined_GC_MD(i, j, k) = sum(sum(BLOCK))/length(a);
     end    
     
     % Excitatory 
     Ln1 = length(l1);
     
     if (Ln1~=0)
       for w = 1:Ln1
          SBLOCK = SBLOCK + BLOCK( l1(w), l2(w) );
       end
     
       EXC_GC_MD (i, j, k) = SBLOCK/Ln1;
     end
     
     
     % Inhibitory 
     SBLOCK = 0;
     
     Ln1 = length(r1);
     
     if (Ln1~=0)
       for w = 1:Ln1
          SBLOCK = SBLOCK + BLOCK( r1(w), r2(w) );
       end
     
       INH_GC_MD (i, j, k) = SBLOCK/Ln1;
     end
     
     
     mj1 = mj1+Nb(j);
    
     
    end   
     
    mi1 = mi1+Nb(i);
    
  end
    
  TEMP0 = Combined_GC_MD(:,:,k); 
  TEMP1 = EXC_GC_MD(:,:,k);
  TEMP2 = INH_GC_MD(:,:,k);   
  
  VC_MD = [VC_MD TEMP0(:)];
  VE_MD = [VE_MD TEMP1(:)];
  VI_MD = [VI_MD TEMP2(:)];
    
    
    
end


GCM_NC = mean(GC_NC, 3);
GCM_MD = mean(GC_MD, 3);

W_EXC_NC = mean(EXC_GC_NC, 3);
W_EXC_MD = mean(EXC_GC_MD, 3);

W_INH_NC = mean(INH_GC_NC, 3);
W_INH_MD = mean(INH_GC_MD, 3);

W_COMBINE_NC = mean(Combined_GC_NC, 3);
W_COMBINE_MD = mean(Combined_GC_MD, 3);

%=============================

% Excitatory
VEM_NC = VE_NC';
VEM_MD = VE_MD';

WEm_NC = mean(VEM_NC);
WEm_MD = mean(VEM_MD);

SEe_NC = std(VEM_NC)/sqrt(n1);
SEe_MD = std(VEM_MD)/sqrt(n2);

WEm_NC_MD = [WEm_NC' WEm_MD'];
SEe_NC_MD = [SEe_NC' SEe_MD'];

% Inhibitory
VIM_NC = VI_NC';
VIM_MD = VI_MD';

WIm_NC = mean(VIM_NC);
WIm_MD = mean(VIM_MD);

SEi_NC = std(VIM_NC)/sqrt(n1);
SEi_MD = std(VIM_MD)/sqrt(n2);

WIm_NC_MD = [WIm_NC' WIm_MD'];
SEi_NC_MD = [SEi_NC' SEi_MD'];


% Combined
VCM_NC = VC_NC';
VCM_MD = VC_MD';

Wm_NC = mean(VCM_NC);
Wm_MD = mean(VCM_MD);

SE_NC = std(VCM_NC)/sqrt(n1);
SE_MD = std(VCM_MD)/sqrt(n2);

Wm_NC_MD = [Wm_NC' Wm_MD'];
SE_NC_MD = [SE_NC' SE_MD'];

%% ========================================================================

if (Flag == 1)
   save GC_NC_MCI.mat GC_NC GC_MD n1 n2
elseif (Flag == 2)
   save GC_NC_AD.mat  GC_NC GC_MD n1 n2
else
   save GC_MCI_AD.mat GC_NC GC_MD n1 n2
end



%% ========================================================================

CLIM1 = [0.5 1.2];
CLIM2 = [-1.2 -0.5];
CLIM3 = [-0.05 0.5];

xvalues = {'DMN', 'SAL', 'EXE', 'LIM'};
yvalues = xvalues;


figure;
h1 = heatmap(xvalues, yvalues, W_COMBINE_NC, 'colormap', parula,'ColorLimits', CLIM3, 'FontSize',14);
h1.CellLabelFormat = '%.2g';
colorbar;
if (Flag==1||Flag==2)
  title('NC');
else
  title('MCI');
end

figure;
h1 = heatmap(xvalues, yvalues, W_COMBINE_MD, 'colormap', parula, 'ColorLimits', CLIM3,'FontSize',14);
if (Flag==1)
  title('MCI');
else
  title('AD');
end
h1.CellLabelFormat = '%.2g';
colorbar;



% figure;
% % h1 = heatmap(xvalues, yvalues, MA_CON_MEAN, 'colormap',  parula, 'CellLabelColor','none', 'FontSize',14 );
% h1 = heatmap(xvalues, yvalues, W_EXC_NC, 'colormap', parula, 'ColorLimits', CLIM1,'FontSize',14);
% title('Excitatory EC (NC)');
% % title('Inhibitory EC (NC)');
% colorbar;
% 
% 
% figure;
% h1 = heatmap(xvalues, yvalues, W_EXC_MD, 'colormap', parula, 'ColorLimits', CLIM1,'FontSize',14);
% title('Excitatory EC (AD)');
% % title('Inhibitory EC (NC)');
% colorbar;
% 
% 
% figure;
% h1 = heatmap(xvalues, yvalues, W_INH_NC, 'colormap', parula, 'ColorLimits', CLIM2,'FontSize',14);
% title('Inhibitory EC (NC)');
% colorbar;
% 
% figure;
% h1 = heatmap(xvalues, yvalues, W_INH_MD, 'colormap', parula, 'ColorLimits', CLIM2,'FontSize',14);
% title('Inhibitory EC (AD)');
% colorbar;
% 
% 
% 
% 







