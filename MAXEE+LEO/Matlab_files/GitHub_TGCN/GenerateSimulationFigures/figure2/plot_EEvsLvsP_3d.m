%This Matlab script can be used to generate Fig. 2 in the article:
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

mm = length(Mrange);

%Compute UL sum SE using the fractions of UL data
sumSE_MMMSE_ULavgss = ULfraction*squeeze(mean(sumSE_MMMSE(:,:,mm,:),4)); %+ DLfraction*sumSE_MMMSE_DL;
sumSE_SMMSE_ULavg = ULfraction*squeeze(mean(sumSE_SMMSE(:,:,mm,:),4)); %+ DLfraction*sumSE_SMMSE_DL;
sumSE_RZF_ULavgss = ULfraction*squeeze(mean(sumSE_RZF(:,:,mm,:),4)); %+ DLfraction*sumSE_RZF_DL;
sumSE_ZF_ULavgss = ULfraction*squeeze(mean(sumSE_ZF(:,:,mm,:),4)); %+ DLfraction*sumSE_ZF_DL;
sumSE_MR_ULavgss = ULfraction*squeeze(mean(sumSE_MR(:,:,mm,:),4)); %+ DLfraction*sumSE_MR_DL;


%% Numerical

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
        [P_MR_UL,P_ZF_UL,P_MMMSE_UL,cp] = functionCPcomputationvsLvsP(lambdaRange,fRange,Mrange,K,L,B,tau_c,tau_p,ULfraction,valueset,sumSE_MR_ULavgss,sumSE_ZF_ULavgss,sumSE_MMMSE_ULavgss);
        
        %Effective transmit power per UE coefficient in W
        Ucalss(l,s) = P0*gamma(alpha/2+1)/Pathlossat1km/mu_UE/(pi*lambdaRange(l))^(alpha/2);
        
        %Compute total effective transmit power
        ETP_totalss(l,s) =  K*Ucalss(l,s)*(tau_p + (tau_c-tau_p)*(ULfraction))/tau_c;
        
        %Compute EE with MR
        EE_MR_ULss(l,s) = (B*sumSE_MR_ULavgss(l,s))./(ETP_totalss(l,s) + P_MR_UL(l,s));
        
        %Compute EE with RZF
        EE_ZF_ULss(l,s) = (B*sumSE_ZF_ULavgss(l,s))./(ETP_totalss(l,s) + P_ZF_UL(l,s));
        
        %Compute EE with M-MMSE
        EE_MMMSE_ULss(l,s) = (B*sumSE_MMMSE_ULavgss(l,s))./(ETP_totalss(l,s) + P_MMMSE_UL(l,s));
        
    end
    
end

% %Plot CP per cell in W vs BS density (\lambda)
% %Circuit power consumption terms
% P_SAME = cp.P_SAME;
% P_CE = cp.P_CE;
% P_CD_MR = cp.P_CD_MR;
% P_BH_MR = cp.P_BH_MR;
% P_SP_UL_MR = cp.P_SP_UL_MR;
% mm = length(Mrange);
% ff = length(fRange);    
% figure(1);
% box on;
% semilogy(lambdaRange,P_SAME(:,ff),'bs-','LineWidth',1); hold on;
% semilogy(lambdaRange,P_CE(:,ff),'r--','LineWidth',1);
% semilogy(lambdaRange,P_CD_MR(:,ff),'k.-','LineWidth',1);
% semilogy(lambdaRange,P_BH_MR(:,ff),'g+-','LineWidth',1);
% semilogy(lambdaRange,P_SP_UL_MR(:,ff),'c-','LineWidth',1);
% semilogy(lambdaRange,ETP_total(:,ff),'md-','LineWidth',1);
% xlabel('BSs density (\lambda)');
% ylabel('CP per cell [W/cell]');
% legend('SAME','ChannEst','CodDec','Backh','SigProc','PowTx','Location','NorthWest');
% title(['CP(\lambda) with \zeta =' num2str(fRange(ff)) ', K = ' num2str(K) '  and M =' num2str(Mrange(mm))])
%     
% 
%% Plot simulation results

[ff, ll] = meshgrid(fRange,lambdaRange);
% interpolation
npointsxi = 10;
npointsyi = 30;
xi = linspace(fRange(1),fRange(end),npointsxi);
yi = linspace(lambdaRange(1),lambdaRange(end),npointsyi);
[fRangeInterpss, lambdaRangeInterpss] = meshgrid(xi, yi);


%Plot Figure 5.14a
figure(1);
hold on; box on; grid on;

EE_MMMSE_ULinterpss = interp2(fRange, lambdaRange, EE_MMMSE_ULss/10^6, fRangeInterpss, lambdaRangeInterpss, 'spline');
surfc(fRangeInterpss, lambdaRangeInterpss, EE_MMMSE_ULinterpss, 'LineStyle', 'none', 'FaceColor', 'interp')
colormap(autumn)
hold on
contour3(fRangeInterpss,lambdaRangeInterpss,EE_MMMSE_ULinterpss,10,'k')
[~,I] = max(EE_MMMSE_ULinterpss(:));
[row,col] = ind2sub(size(EE_MMMSE_ULinterpss),I); % 2D maximizer
fopt = fRangeInterpss(1,col);  % Nms maximizer(Users=fixed)
lopt = lambdaRangeInterpss(row,1);  % Nbs maximizer (Users=fixed)
hold on
plot3(fopt,lopt,EE_MMMSE_ULinterpss(row,col),'k*','MarkerSize',16,'MarkerFaceColor','black');
hold on
plot3(fopt,lopt,min(min(EE_MMMSE_ULinterpss)),'k*','MarkerSize',16,'MarkerFaceColor','black');
zlim([min(min(EE_MMMSE_ULinterpss)) 9])

view([-17 32]);

% xlabel('Pilot reuse factor (\zeta)')
% ylabel('BS density (\lambda)');
% zlabel('EE [Mbit/Joule]');
% title('EE(\lambda,\zeta) with M-MMSE combiner')
xlabel('zeta')
ylabel('lambda');
zlabel('EE');
set(gca,'FontSize',20)


%Plot Figure 5.14b
figure(2);
hold on; box on; grid on;

EE_ZF_ULinterpss = interp2(fRange, lambdaRange, EE_ZF_ULss/10^6, fRangeInterpss, lambdaRangeInterpss, 'spline');
surfc(fRangeInterpss, lambdaRangeInterpss, EE_ZF_ULinterpss, 'LineStyle', 'none', 'FaceColor', 'interp')
colormap(autumn)
hold on
contour3(fRangeInterpss,lambdaRangeInterpss,EE_ZF_ULinterpss,10,'k')
[~,I] = max(EE_ZF_ULinterpss(:));
[row,col] = ind2sub(size(EE_ZF_ULinterpss),I); % 2D maximizer
fopt = fRangeInterpss(1,col);  % Nms maximizer(Users=fixed)
lopt = lambdaRangeInterpss(row,1);  % Nbs maximizer (Users=fixed)
hold on
plot3(fopt,lopt,EE_ZF_ULinterpss(row,col),'k*','MarkerSize',16,'MarkerFaceColor','black');
hold on
plot3(fopt,lopt,min(min(EE_ZF_ULinterpss)),'k*','MarkerSize',16,'MarkerFaceColor','black');
zlim([min(min(EE_ZF_ULinterpss)) 9])

view([-17 32]);

% xlabel('Pilot reuse factor (\zeta)')
% ylabel('BS density (\lambda)');
% zlabel('EE [Mbit/Joule]');
% title('EE(\lambda,\zeta) with ZF combiner')
xlabel('zeta')
ylabel('lambda');
zlabel('EE');
set(gca,'FontSize',20)


%Plot Figure 5.14c
figure(3);
hold on; box on; grid on;

EE_MR_ULinterpss = interp2(fRange, lambdaRange, EE_MR_ULss/10^6, fRangeInterpss, lambdaRangeInterpss, 'spline');
surfc(fRangeInterpss, lambdaRangeInterpss, EE_MR_ULinterpss, 'LineStyle', 'none', 'FaceColor', 'interp')
colormap(autumn)
hold on
contour3(fRangeInterpss,lambdaRangeInterpss,EE_MR_ULinterpss,10,'k')
[~,I] = max(EE_MR_ULinterpss(:));
[row,col] = ind2sub(size(EE_MR_ULinterpss),I); % 2D maximizer
fopt = fRangeInterpss(1,col);  % Nms maximizer(Users=fixed)
lopt = lambdaRangeInterpss(row,1);  % Nbs maximizer (Users=fixed)
hold on
plot3(fopt,lopt,EE_MR_ULinterpss(row,col),'k*','MarkerSize',16,'MarkerFaceColor','black');
hold on
plot3(fopt,lopt,min(min(EE_MR_ULinterpss)),'k*','MarkerSize',16,'MarkerFaceColor','black');
zlim([min(min(EE_MR_ULinterpss)) 9])

view([-17 32]);

% xlabel('Pilot reuse factor (\zeta)')
% ylabel('BS density (\lambda)');
% zlabel('EE [Mbit/Joule]');
% title('EE(\lambda,\zeta) with MR combiner')
xlabel('zeta')
ylabel('lambda');
zlabel('EE');
set(gca,'FontSize',20)


%% Numerical Multi Slope

%Load SE simulation data related to the multislope path loss model that are 
%generated using the code in ComputeULAvgErgodicSE
load ../../SimulationResults/sumSEULvsLPwithPA_3slopes_lambda1to60.mat

% M = 100
mm = length(Mrange);

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

% %Plot CP per cell in W vs BS density (\lambda)
% %Circuit power consumption terms
% P_SAME = cp.P_SAME;
% P_CE = cp.P_CE;
% P_CD_MR = cp.P_CD_MR;
% P_BH_MR = cp.P_BH_MR;
% P_SP_UL_MR = cp.P_SP_UL_MR;
% mm = length(Mrange);
% ff = length(fRange);    
% figure(1);
% box on;
% semilogy(lambdaRange,P_SAME(:,ff),'bs-','LineWidth',1); hold on;
% semilogy(lambdaRange,P_CE(:,ff),'r--','LineWidth',1);
% semilogy(lambdaRange,P_CD_MR(:,ff),'k.-','LineWidth',1);
% semilogy(lambdaRange,P_BH_MR(:,ff),'g+-','LineWidth',1);
% semilogy(lambdaRange,P_SP_UL_MR(:,ff),'c-','LineWidth',1);
% semilogy(lambdaRange,ETP_total(:,ff),'md-','LineWidth',1);
% 
% xlabel('BSs density (\lambda)');
% ylabel('CP per cell [W/cell]');
% legend('SAME','ChannEst','CodDec','Backh','SigProc','PowTx','Location','NorthWest');
% title(['CP(\lambda) with \zeta =' num2str(fRange(ff)) ', K = ' num2str(K) '  and M =' num2str(Mrange(mm))])
% % ylim([0 16]);
    

%% Plot simulation results

[ff, ll] = meshgrid(fRange,lambdaRange);
% interpolation
npointsxi = 10;
npointsyi = 30;
xi = linspace(fRange(1),fRange(end),npointsxi);
yi = linspace(lambdaRange(1),lambdaRange(end),npointsyi);
[fRangeInterpms, lambdaRangeInterpms] = meshgrid(xi, yi);


%Plot Figure 5.14a
figure(11);
hold on; box on; grid on;

EE_MMMSE_ULinterpms = interp2(fRange, lambdaRange, EE_MMMSE_ULms/10^6, fRangeInterpms, lambdaRangeInterpms, 'spline');
surfc(fRangeInterpms, lambdaRangeInterpms, EE_MMMSE_ULinterpms, 'LineStyle', 'none', 'FaceColor', 'interp')
colormap(autumn)
hold on
contour3(fRangeInterpms,lambdaRangeInterpms,EE_MMMSE_ULinterpms,10,'k')
[~,I] = max(EE_MMMSE_ULinterpms(:));
[row,col] = ind2sub(size(EE_MMMSE_ULinterpms),I); % 2D maximizer
fopt = fRangeInterpms(1,col);  % Nms maximizer(Users=fixed)
lopt = lambdaRangeInterpms(row,1);  % Nbs maximizer (Users=fixed)
hold on
plot3(fopt,lopt,EE_MMMSE_ULinterpms(row,col),'k*','MarkerSize',16,'MarkerFaceColor','black');
hold on
plot3(fopt,lopt,min(min(EE_MMMSE_ULinterpms)),'k*','MarkerSize',16,'MarkerFaceColor','black');
zlim([min(min(EE_MMMSE_ULinterpms)) 9])

view([245 32]);
% xlabel('Pilot reuse factor (\zeta)')
% ylabel('BS density (\lambda)');
% zlabel('EE [Mbit/Joule]');
% title('EE(\lambda,\zeta) with M-MMSE combiner')
xlabel('zeta')
ylabel('lambda');
zlabel('EE');
set(gca,'FontSize',20)


%Plot Figure 5.14b
figure(22);
hold on; box on; grid on;

EE_ZF_ULinterpms = interp2(fRange, lambdaRange, EE_ZF_ULms/10^6, fRangeInterpms, lambdaRangeInterpms, 'spline');
surfc(fRangeInterpms, lambdaRangeInterpms, EE_ZF_ULinterpms, 'LineStyle', 'none', 'FaceColor', 'interp')
colormap(autumn)
hold on
contour3(fRangeInterpms,lambdaRangeInterpms,EE_ZF_ULinterpms,10,'k')
[~,I] = max(EE_ZF_ULinterpms(:));
[row,col] = ind2sub(size(EE_ZF_ULinterpms),I); % 2D maximizer
fopt = fRangeInterpms(1,col);  % Nms maximizer(Users=fixed)
lopt = lambdaRangeInterpms(row,1);  % Nbs maximizer (Users=fixed)
hold on
plot3(fopt,lopt,EE_ZF_ULinterpms(row,col),'k*','MarkerSize',16,'MarkerFaceColor','black');
hold on
plot3(fopt,lopt,min(min(EE_ZF_ULinterpms)),'k*','MarkerSize',16,'MarkerFaceColor','black');
zlim([min(min(EE_ZF_ULinterpms)) 9])

view([245 32]);

% xlabel('Pilot reuse factor (\zeta)')
% ylabel('BS density (\lambda)');
% zlabel('EE [Mbit/Joule]');
% title('EE(\lambda,\zeta) with ZF combiner')
xlabel('zeta')
ylabel('lambda');
zlabel('EE');
set(gca,'FontSize',20)


%Plot Figure 5.14c
figure(33);
hold on; box on; grid on;

EE_MR_ULinterpms = interp2(fRange, lambdaRange, EE_MR_ULms/10^6, fRangeInterpms, lambdaRangeInterpms, 'spline');
surfc(fRangeInterpms, lambdaRangeInterpms, EE_MR_ULinterpms, 'LineStyle', 'none', 'FaceColor', 'interp')
colormap(autumn)
hold on
contour3(fRangeInterpms,lambdaRangeInterpms,EE_MR_ULinterpms,10,'k')
[~,I] = max(EE_MR_ULinterpms(:));
[row,col] = ind2sub(size(EE_MR_ULinterpms),I); % 2D maximizer
fopt = fRangeInterpms(1,col);  % Nms maximizer(Users=fixed)
lopt = lambdaRangeInterpms(row,1);  % Nbs maximizer (Users=fixed)
hold on
plot3(fopt,lopt,EE_MR_ULinterpms(row,col),'k*','MarkerSize',16,'MarkerFaceColor','black');
hold on
plot3(fopt,lopt,min(min(EE_MR_ULinterpms)),'k*','MarkerSize',16,'MarkerFaceColor','black');
zlim([min(min(EE_MR_ULinterpms)) 9])

view([245 32]);
% xlabel('Pilot reuse factor (\zeta)')
% ylabel('BS density (\lambda)');
% zlabel('EE [Mbit/Joule]');
% title('EE(\lambda,\zeta) with MR combiner')
xlabel('zeta')
ylabel('lambda');
zlabel('EE');
set(gca,'FontSize',20)