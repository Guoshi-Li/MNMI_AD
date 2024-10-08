
% Function to compute the BOLD signals 
% Need to call a SPM function ("spm_hrf.m")

function [BOLD] = Model_HEMO_HRF(N, dt, ST, XX, FLAG_Mean_BOLD) 

NR = N; % Number of brain regions
DT = dt;   % Integration step (sec)
NE = 4;    % Number of equations

%=================================== 
%     Hemodynamic model
%=================================== 
  
T = 16;
TS = ST;
p = [6 16 1 1 6 0 TS];
RT = DT; 
[L1, L2] = size(XX);

% Create the HRF by calling the SPM function
hrf = spm_hrf(RT, p, T); 

 for i = 1:NR
    
   E = XX(:,i);
   I = XX(:,i+NR);
   x = (2/3)*E + (1/3)*I;
  
   if (FLAG_Mean_BOLD==1)
     XM = mean(x);
     x = x - XM;
   end
  
   % Compute the convolution of composite neural activity and HRF 
   y = conv(x, hrf);
  
   BOLD(:,i) = y(1:L1);
  

 end

end



