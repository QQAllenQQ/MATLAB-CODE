%This Matlab script can be used to generate Fig. 4 in the article:
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

%M = 100;
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
[~,col_MRss] = ind2sub(size(EE_MR_ULss),I); % 2D maximizer
fopt_MRss = fRange(1,col_MRss);  % Nms maximizer(Users=fixed)
%ZF
[~,I] = max(EE_ZF_ULss(:));
[~,col_ZFss] = ind2sub(size(EE_ZF_ULss),I); % 2D maximizer
fopt_ZFss = fRange(1,col_ZFss);  % Nms maximizer(Users=fixed)
%M-MMSE
[~,I] = max(EE_MMMSE_ULss(:));
[~,col_MMMSEss] = ind2sub(size(EE_MMMSE_ULss),I); % 2D maximizer
fopt_MMMSEss = fRange(1,col_MMMSEss);  % Nms maximizer(Users=fixed)

% Single-slope EE vs BS density
figure(1);
hold on; box on;
mmmse_ss = plot(lambdaRange,smooth(EE_MMMSE_ULss(:,fopt_MMMSEss)/1e6,'sgolay',3),'rd','LineWidth',2,'MarkerSize',12);
zf_ss = plot(lambdaRange,smooth(EE_ZF_ULss(:,fopt_ZFss)/1e6,'sgolay',3),'ko','LineWidth',2,'MarkerSize',12);
mr_ss = plot(lambdaRange,smooth(EE_MR_ULss(:,fopt_MRss)/1e6,'sgolay',3),'bs','LineWidth',2,'MarkerSize',12);
plot(lambdaRange,smooth(EE_ZF_ULss(:,fopt_ZFss)/1e6,'sgolay',3),'ko-','LineWidth',2,'MarkerSize',12);
plot(lambdaRange,smooth(EE_MR_ULss(:,fopt_MRss)/1e6,'sgolay',3),'bs-','LineWidth',2,'MarkerSize',12);
mmmse2_ss = plot(lambdaRange,smooth(EE_MMMSE_ULss(:,fopt_MMMSEss)/1e6,'sgolay',3),'r-','LineWidth',2,'MarkerSize',12);
legend('M-MMSE', 'ZF', 'MR', 'Location', 'NorthEast')
xlabel('BS density (\lambda) [BSs/km^2]')
ylabel('EE [bit/Joule]')
title('(a) EE versus BS density (\lambda)')

set(gca, 'FontSize', 14)
grid on;
xticks(0:5:max(lambdaRange))

%% Numerical Multi Slope

%Load SE simulation data related to the mulitislope path loss model that are 
%generated using the code in ComputeULAvgErgodicSE
load ../../SimulationResults/sumSEULvsLPwithPA_3slopes_lambda1to60.mat

%Compute UL sum SE using the fractions of UL data
sumSE_MMMSE_ULavgms = ULfraction*squeeze(mean(sumSE_MMMSE(:,:,mm,:),4)); %+ DLfraction*sumSE_MMMSE_DL;
sumSE_SMMSE_ULavgms = ULfraction*squeeze(mean(sumSE_SMMSE(:,:,mm,:),4)); %+ DLfraction*sumSE_SMMSE_DL;
sumSE_RZF_ULavgms = ULfraction*squeeze(mean(sumSE_RZF(:,:,mm,:),4)); %+ DLfraction*sumSE_RZF_DL;
sumSE_ZF_ULavgms = ULfraction*squeeze(mean(sumSE_ZF(:,:,mm,:),4)); %+ DLfraction*sumSE_ZF_DL;
sumSE_MR_ULavgms = ULfraction*squeeze(mean(sumSE_MR(:,:,mm,:),4)); %+ DLfraction*sumSE_MR_DL;


%Prepare to save simulation results
EE_MR_ULms = zeros(length(lambdaRange),length(fRange));
EE_ZF_ULms = zeros(length(lambdaRange),length(fRange));
EE_MMMSE_ULms = zeros(length(lambdaRange),length(fRange));
ETP_totalms = zeros(length(lambdaRange),length(fRange));
Ucalms = zeros(length(lambdaRange),length(fRange));

%Multi-slope model Pathloss exponent
pathLossExp_vec = Pathloss_model.pathLossExp_vec;
distancesPathLossExp_vec = Pathloss_model.distancesPathLossExp_vec;
Gn_vec = Pathloss_model.Gn_vec;

% Go through all number of BS densities
for l = 1:length(lambdaRange)
    
    lambda = lambdaRange(l);
    % Number of BSs
    L = lambda*(squareLength/1000)^2;
    
    %Go through all number of Pilot reuse factors
    for s = 1:length(fRange)
        
        %Compute length of pilot sequences
        tau_p = K*fRange(s);
        
        %Compute the total CP with different schemes
        [P_MR_UL,P_ZF_UL,P_MMMSE_UL,cp] = functionCPcomputationvsLvsP(lambdaRange,fRange,Mrange,K,L,B,tau_c,tau_p,ULfraction,valueset,sumSE_MR_ULavgms,sumSE_ZF_ULavgms,sumSE_MMMSE_ULavgms);
        
        %Effective transmit power per UE coefficient in W
        Ucalms_sum = 0;
        for aa = 1:length(pathLossExp_vec)
            alphan = pathLossExp_vec(end-aa+1);
            dn = distancesPathLossExp_vec(aa)/1e3;
            dnplus1 = distancesPathLossExp_vec(aa+1)/1e3;
            gamma_difference = ( gammainc(pi*lambda*dn^2,1+alphan/2,'upper')-gammainc(pi*lambda*dnplus1^2,1+alphan/2,'upper') )*gamma(1+alphan/2);
            Ucal_term = gamma_difference/Gn_vec(end-aa+1)/(pi*lambdaRange(l))^(alphan/2);
            Ucalms_sum = Ucalms_sum + Ucal_term;
        end
        Ucalms(l,s) = P0/mu_UE*Ucalms_sum;
        
        %         Ucalss(l,s) = gamma(alpha/2+1)*P0/mu_UE/Pathlossat1km/(pi*lambdaRange(l))^(alpha/2);
        %Compute total effective transmit power
        ETP_totalms(l,s) =  K*Ucalms(l,s)*(tau_p + (tau_c-tau_p)*(ULfraction))/tau_c;
        
        %Compute EE with MR
        EE_MR_ULms(l,s) = (B*sumSE_MR_ULavgms(l,s))./(ETP_totalms(l,s) + P_MR_UL(l,s));
        
        %Compute EE with RZF
        EE_ZF_ULms(l,s) = (B*sumSE_ZF_ULavgms(l,s))./(ETP_totalms(l,s) + P_ZF_UL(l,s));
        
        %Compute EE with M-MMSE
        EE_MMMSE_ULms(l,s) = (B*sumSE_MMMSE_ULavgms(l,s))./(ETP_totalms(l,s) + P_MMMSE_UL(l,s));
        
    end
    
end


%% Compute the optimal pilot reuse factor
% MR
[~,I] = max(EE_MR_ULms(:));
[~,col_MRms] = ind2sub(size(EE_MR_ULms),I); % 2D maximizer
fopt_MRms = fRange(1,col_MRms);  % Nms maximizer(Users=fixed)
%ZF
[~,I] = max(EE_ZF_ULms(:));
[~,col_ZFms] = ind2sub(size(EE_ZF_ULms),I); % 2D maximizer
fopt_ZFms = fRange(1,col_ZFms);  % Nms maximizer(Users=fixed)
%M-MMSE
[~,I] = max(EE_MMMSE_ULms(:));
[~,col_MMMSEms] = ind2sub(size(EE_MMMSE_ULms),I); % 2D maximizer
fopt_MMMSEms = fRange(1,col_MMMSEms);  % Nms maximizer(Users=fixed)

%%
% Multi-slope EE vs BS density
mmmse_ms = plot(lambdaRange,smooth(EE_MMMSE_ULms(:,col_MMMSEms)/1e6,'sgolay',3),'r--','LineWidth',2,'MarkerSize',12);
plot(lambdaRange,smooth(EE_MMMSE_ULms(:,col_MMMSEms)/1e6,'sgolay',3),'rd--','LineWidth',2,'MarkerSize',12);
plot(lambdaRange,smooth(EE_ZF_ULms(:,col_ZFms)/1e6,'sgolay',3),'ko--','LineWidth',2,'MarkerSize',12);
plot(lambdaRange,smooth(EE_MR_ULms(:,col_MRms)/1e6,'sgolay',3),'bs--','LineWidth',2,'MarkerSize',12);

legend('M-MMSE', 'ZF', 'MR', 'Location', 'NorthEast')
xlabel('BS density (\lambda) [BSs/km^2]')
ylabel('EE [bit/Joule]')
title('(a) EE versus BS density (\lambda)')

set(gca, 'FontSize', 14)
grid on;
xticks(0:5:max(lambdaRange))
