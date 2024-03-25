%This Matlab script can be used to generate the single slope version of the
%Fig 4 (multislope) in the article:
%
% Andrea Pizzo, Daniel Verenzuela, Luca Sanguinetti and Emil Björnson, "Network Deployment for Maximal Energy Efficiency 
% in Uplink with Multislope Path Loss," IEEE Transactions on Green Communications and Networking, Submitted to.
%
%This is version 1.0 (Last edited: 10-April-2018)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%

%Initialize
close all;
clear;
clc;

%Select the value set of the power coefficients that is considered 
%(set valueset=2 to use the hardware parameters listen in Table 1 and 
% generate the results as in the article)
valueset = 2;

%% Numerical Single Slope

%Load SE simulation data related to the 1 slope path loss model that are 
%generated using the code in ComputeULAvgErgodicSE
load ../../SimulationResults/sumSEULvsLPwithPA_1slope_lambda1to60.mat

%Load fixed Propagation parameters
run ../../ComputeULAvgErgodicSE/SetPropagationParameters.m;

%Retrieve the SNRs value in dB in the inverse statistical power control policy
SNR = 10^(DeltadB/10);
SNRp = 10^(DeltadB_pilot/10);

% M = 100
mm = length(Mrange);

%Compute UL sum SE using the fractions of UL data
sumSE_MMMSE_ULavg = ULfraction*squeeze(mean(sumSE_MMMSE(:,:,mm,:),4)); %+ DLfraction*sumSE_MMMSE_DL;
sumSE_SMMSE_ULavg = ULfraction*squeeze(mean(sumSE_SMMSE(:,:,mm,:),4)); %+ DLfraction*sumSE_SMMSE_DL;
sumSE_RZF_ULavg = ULfraction*squeeze(mean(sumSE_RZF(:,:,mm,:),4)); %+ DLfraction*sumSE_RZF_DL;
sumSE_ZF_ULavg = ULfraction*squeeze(mean(sumSE_ZF(:,:,mm,:),4)); %+ DLfraction*sumSE_ZF_DL;
sumSE_MR_ULavg = ULfraction*squeeze(mean(sumSE_MR(:,:,mm,:),4)); %+ DLfraction*sumSE_MR_DL;

%Prepare to save simulation results
EE_MR_ULss = zeros(length(lambdaRange),length(fRange));
EE_ZF_ULss = zeros(length(lambdaRange),length(fRange));
EE_MMMSE_ULss = zeros(length(lambdaRange),length(fRange));
ETP_totalss = zeros(length(lambdaRange),length(fRange));
Ucalss = zeros(length(lambdaRange),length(fRange));

% Go through all number of BS densities
for l = 1:length(lambdaRange)
    
    % Number of BSs
    L = lambdaRange(l)*(squareLength/1000)^2;
    
    %Go through all number of Pilot reuse factors
    for s = 1:length(fRange)
        
        %Compute length of pilot sequences
        tau_p = K*fRange(s);
        
        %Compute the total CP with different schemes
        [P_MR_UL,P_ZF_UL,P_MMMSE_UL,cp] = functionCPcomputationvsLvsP(lambdaRange,fRange,Mrange,K,L,B,tau_c,tau_p,ULfraction,valueset,sumSE_MR_ULavg,sumSE_ZF_ULavg,sumSE_MMMSE_ULavg);
        
        %Effective transmit power per UE coefficient in W
        Ucalss(l,s) = P0*gamma(alpha/2+1)/Pathlossat1km/mu_UE/(pi*lambdaRange(l))^(alpha/2);
        
        %Compute total effective transmit power
        ETP_totalss(l,s) =  K*Ucalss(l,s)*(tau_p + (tau_c-tau_p)*(ULfraction))/tau_c;
        
        
        %Compute EE with MR
        EE_MR_ULss(l,s) = (B*sumSE_MR_ULavg(l,s))./(ETP_totalss(l,s) + P_MR_UL(l,s));
        
        %Compute EE with RZF
        EE_ZF_ULss(l,s) = (B*sumSE_ZF_ULavg(l,s))./(ETP_totalss(l,s) + P_ZF_UL(l,s));
        
        %Compute EE with M-MMSE
        EE_MMMSE_ULss(l,s) = (B*sumSE_MMMSE_ULavg(l,s))./(ETP_totalss(l,s) + P_MMMSE_UL(l,s));
        
    end
    
end

%% Compute the optimal pilot reuse factor
% MR
[~,I] = max(EE_MR_ULss(:));
[~,col_MRss_num] = ind2sub(size(EE_MR_ULss),I); % 2D maximizer
fopt_MRss_num = fRange(1,col_MRss_num);  % Nms maximizer(Users=fixed)
%ZF
[~,I] = max(EE_ZF_ULss(:));
[~,col_ZFss_num] = ind2sub(size(EE_ZF_ULss),I); % 2D maximizer
fopt_ZFss_num = fRange(1,col_ZFss_num);  % Nms maximizer(Users=fixed)
%M-MMSE
[~,I] = max(EE_MMMSE_ULss(:));
[~,col_MMMSEss_num] = ind2sub(size(EE_MMMSE_ULss),I); % 2D maximizer
fopt_MMMSEss_num = fRange(1,col_MMMSEss_num);  % Nms maximizer(Users=fixed)


%% Theoretical (single-slope path-loss model)

%Compute theoretical aggregate SE per cell vs BS density vs pilot reuse factor (only MR and ZF)
sumSE_MR_UL_theory = NaN*ones(length(lambdaRange),length(fRange),length(Mrange));
sumSE_ZF_UL_theory = NaN*ones(length(lambdaRange),length(fRange),length(Mrange));

% MR
%Go through all the BS densities
for l = 1:length(lambdaRange)
    lambda = lambdaRange(l);
    %Go through all pilot reuse factors
    for s = 1:length(fRange)
        f = fRange(s);
        tau_p = K*f;
        %Go through all number of antennas
        for m = 1:length(Mrange)
            M = Mrange(m);
            interf_PC = M/(alpha-1)/f;
            interf_noPC_term1 = (K+1/SNR)*(1+2/(alpha-2)/f+1/SNRp);
            interf_noPC_term2 = 2*K/(alpha-2)*(1+1/SNRp);
            interf_noPC_term3 = K/f*(4/(alpha-2)^2+1/(alpha-1));
            interf_noPC = interf_noPC_term1 + interf_noPC_term2 + interf_noPC_term3;
            SINR_theory = M/(interf_PC + interf_noPC);
            %save
            sumSE_MR_UL_theory(l,s,m) = ULfraction*K*(1-tau_p/tau_c)*log2(1 + SINR_theory);
        end
    end
end
% ZF
%Go through all the BS densities
for l = 1:length(lambdaRange)
    lambda = lambdaRange(l);
    %Go through all pilot reuse factor
    for s = 1:length(fRange)
        f = fRange(s);
        tau_p = K*f;
        %Go through all number of antennas
        for m = 1:length(Mrange)
            M = Mrange(m);
            interf_PC = (M-K)/(alpha-1)/f;
            interf_noPC_term1 = (K+1/SNR)*(1+2/(alpha-2)/f+1/SNRp);
            interf_noPC_term2 = 2*K/(alpha-2)*(1+1/SNRp);
            interf_noPC_term3 = K/f*(4/(alpha-2)^2+1/(alpha-1));
            interf_noPC_term4 = K*(1+1/(alpha-1)/f);
            interf_noPC = interf_noPC_term1 + interf_noPC_term2 + interf_noPC_term3 - interf_noPC_term4;
            SINR_theory = (M-K)/(interf_PC + interf_noPC);
            %save
            sumSE_ZF_UL_theory(l,s,m) = ULfraction*K*(1-tau_p/tau_c)*log2(1 + SINR_theory);
        end
    end
end

%Prepare to save simulation results
EE_MR_ULss_th = zeros(length(lambdaRange),length(fRange));
EE_ZF_ULss_th = zeros(length(lambdaRange),length(fRange));
ETP_totalss_th = zeros(length(lambdaRange),length(fRange));
Ucalss_th = zeros(length(lambdaRange),length(fRange));

% Go through all number of BS densities
for l = 1:length(lambdaRange)
    
    % Number of BSs
    L = lambdaRange(l)*(squareLength/1000)^2;
    
    %Go through all number of Pilot reuse factors
    for s = 1:length(fRange)
        
        %Compute length of pilot sequences
        tau_p = K*fRange(s);
        
        %Compute the total CP with different schemes
        [P_MR_UL,P_ZF_UL] = functionCPcomputationvsLvsP(lambdaRange,fRange,Mrange,K,L,B,tau_c,tau_p,ULfraction,valueset,sumSE_MR_UL_theory,sumSE_ZF_UL_theory);
        
        %Effective transmit power per UE coefficient in W
        Ucalss(l,s) = P0*gamma(alpha/2+1)/Pathlossat1km/mu_UE/(pi*lambdaRange(l))^(alpha/2);
        
        %Compute total effective transmit power
        ETP_totalss(l,s) =  K*Ucalss(l,s)*(tau_p + (tau_c-tau_p)*(ULfraction))/tau_c;
        
        
        %Compute EE with MR
        EE_MR_ULss_th(l,s) = (B*sumSE_MR_UL_theory(l,s))./(ETP_totalss(l,s) + P_MR_UL(l,s));
        
        %Compute EE with RZF
        EE_ZF_ULss_th(l,s) = (B*sumSE_ZF_UL_theory(l,s))./(ETP_totalss(l,s) + P_ZF_UL(l,s));
        
    end
    
end

% Compute the optimal pilot reuse factor
% MR
[~,I] = max(EE_MR_ULss_th(:));
[~,col_MRms_th] = ind2sub(size(EE_MR_ULss_th),I); % 2D maximizer
fopt_MRms_th = fRange(1,col_MRms_th);  % Nms maximizer(Users=fixed)
%ZF
[~,I] = max(EE_ZF_ULss_th(:));
[~,col_ZFms_th] = ind2sub(size(EE_ZF_ULss_th),I); % 2D maximizer
fopt_ZFms_th = fRange(1,col_ZFms_th);  % Nms maximizer(Users=fixed)

%Plot aggregate SE per cell vs BS density
figure(1);
hold on; box on;
plot(lambdaRange,smooth(EE_MMMSE_ULss(:,fopt_MMMSEss_num)/1e6,'sgolay',2),'rd-','LineWidth',2,'MarkerSize',10);
plot(lambdaRange,smooth(EE_ZF_ULss(:,fopt_ZFss_num)/1e6,'sgolay',2),'ko-','LineWidth',2,'MarkerSize',10);
plot(lambdaRange,smooth(EE_MR_ULss(:,fopt_MRss_num)/1e6,'sgolay',2),'bs-','LineWidth',2,'MarkerSize',10);

%Plot theoretical aggregate SE per cell vs BS density (only MR and ZF)    
figure(1);
hold on; box on;
plot(lambdaRange,smooth(EE_ZF_ULss_th(:,fopt_ZFms_th)/1e6,'sgolay',2),'ko-','LineWidth',2,'MarkerFaceColor',[.5,.5,.5],'MarkerSize',10);
plot(lambdaRange,smooth(EE_MR_ULss_th(:,fopt_MRms_th)/1e6,'sgolay',2),'bo-','LineWidth',2,'MarkerFaceColor',[.5,.5,.5],'MarkerSize',10);

% xlabel('BSs density (\lambda)');
% ylabel('Average sum SE [bit/s/Hz/cell]');
% legend('M-MMSE','S-MMSE','RZF','ZF','MR','ZF-theory','MR-theory','Location','NorthWest');
% title(['sum ASE(\lambda) with \zeta =' num2str(fRange(ff)) ' and M =' num2str(Mrange(mm))])
xlabel('BS density');
ylabel('EE');
legend('MMMSE','ZF','MR','ZFtheory','MRtheory','Location','Best');
grid on;
set(gca,'FontSize',20)
