%This Matlab script can be used to generate the Fig 5 in the article:
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
%This function generate the Fig 5 as it appears in the article. Both ZF and
%MR combining schemes can be used to compute the system energy efficiency.
%

%Initialize
close all;
clear;
clc;

%Select the value set of the power coefficients that is considered 
%(set valueset=2 to use the hardware parameters listen in Table 1 and 
% generate the results as in the article)
valueset = 2;

%Load fixed Propagation parameters
run ../../ComputeULAvgErgodicSE/SetPropagationParameters.m;

%Retrieve the SNRs value in dB in the inverse statistical power control policy
SNR = 10^(DeltadB/10);
SNRp = 10^(DeltadB_pilot/10);

% Generates the large scale path-loss coefficient using (1) single slope or
% (2) multi slope model
pathloss = 2;

%Select the combining scheme at BS using (1) MR or (2) ZF
flagcombScheme = 2;

% Select the optimization algorithm for retrieving the optimal K using (1)
% exhaustive search as in Lemma 4 or (2) simplified as in Corollary 4 in
% the article
flagOptK = 2;

%Define the range of served UEs per cell
Krange = 1:25;

%Define the range of BS antennas
Mrange = 1:250;

%Load simulation settings
define;

%Placeholders for storing of simulation results
ThroughputPerCell = zeros(length(Mrange),length(Krange)); %EE for different lambda and gamma values (using theoretical formulas)
ZetaOpt = zeros(length(Mrange),length(Krange)); %Store optimal beta for each point using the theoretical formulas
PC_total = zeros(length(Mrange),length(Krange));
EE_UL = zeros(length(Mrange),length(Krange));

%Compute the optimal pilot reuse factor
if pathloss == 2    %multislope
    %Compute the auxiliary \mu's terms
    [mu_kappa1, mu_kappa2] = functionComputeDistanceExpectation(pathLossExp_vec,distancesPathLossExp_vec,Gn_vec,lambda);
end

%     %Check for feasibility
%     if gammaval/(alpha-1)>=1
%         error(0,'The achieved rate can be at most: r = log2(alpha)')
%     end

%% THEORETICAL EE

%Go through all number of users
for k = 1:length(Krange)
    
    %select the current number of UEs K
    K = Krange(k);
    
    %Go through all number of antennas
    for m = K:length(Mrange)
        
        %select the current number of BS antennas M
        M = Mrange(m);
        
        %Compute the optimal pilot reuse factor
        if pathloss == 1    %single-slope
            
            %Compute the constants B1 and B2 in Lemma 3
            B1 = K*( 2/(alpha-2) + 4/(alpha-2)^2 - 1/(alpha-1) ) + M/(alpha-1) + 2/(alpha-2)/SNR;
            B2 = K*( 1/SNRp + 2/(alpha-2)*(1 + 1/SNRp) ) + 1/SNR*(1 + 1/SNRp);
            
            if flagcombScheme ==1   %MR combining
                %Compute the optimal pilot reuse factor in Lemma3
                zetaval = (B1*gammaval+2*K/(alpha-1)*gammaval)/(M-K*gammaval-B2*gammaval);
            elseif flagcombScheme == 2  %ZF combining
                %Compute the optimal pilot reuse factor in Lemma3
                zetaval = B1*gammaval/(M-K-B2*gammaval);
            end
            
        elseif pathloss == 2    %multislope
            
            %Compute the constants B1 and B2 in Lemma 3
            B1 = K*( mu_kappa1 + mu_kappa1^2 - mu_kappa2 ) + M*mu_kappa2 + mu_kappa1/SNR;
            B2 = K*( 1/SNRp + mu_kappa1*(1 + 1/SNRp) ) + 1/SNR*(1 + 1/SNRp);
            
            if flagcombScheme ==1   %MR combining
                %Compute the optimal pilot reuse factor in Lemma3
                zetaval = (B1*gammaval+2*K*mu_kappa2*gammaval)/(M-K*gammaval-B2*gammaval);
            elseif flagcombScheme == 2  %ZF combining
                %Compute the optimal pilot reuse factor in Lemma3
                zetaval = B1*gammaval/(M-K-B2*gammaval);
            end
            
        end
        %Save
        ZetaOpt(m,k) = zetaval;
        
        %Check if the two constraints in Eq. (23) are satisfied
        if zetaval>=1 && K*zetaval<=tau_c
            
            %Compute the Area Rate
            ThroughputPerCell(m,k) = B*K*ULfraction*(1-zetaval*K/tau_c) * log2(1+ gammaval);
            
            %Compute the Area Power Consumption APC
            %Compute length of pilot sequences
            tau_p = K*zetaval;
            
            %Compute the total CP with ZF combining scheme
            if flagcombScheme ==1   %MR combining
                [P_CP,cp] = functionCPcomputationvsMvsK_MR(M,K,B,tau_c,tau_p,ULfraction,valueset,ThroughputPerCell(M,K));
            elseif flagcombScheme == 2  %ZF combining
                [P_CP,cp] = functionCPcomputationvsMvsK_ZF(M,K,B,tau_c,tau_p,ULfraction,valueset,ThroughputPerCell(M,K));
            end
            
            %Compute total effective transmit power
            if pathloss == 1    %single-slope
                
                Ucal = P0*gamma(alpha/2+1)/Pathlossat1km/mu_UE/(pi*lambda)^(alpha/2);
                
            elseif pathloss == 2    %multislope
                
                Ucal_sum = 0;
                for aa = 1:length(pathLossExp_vec)
                    alphan = pathLossExp_vec(end-aa+1);
                    dn = distancesPathLossExp_vec(aa)/1e3;
                    dnplus1 = distancesPathLossExp_vec(aa+1)/1e3;
                    gamma_difference = ( gammainc(pi*lambda*dn^2,1+alphan/2,'upper')-gammainc(pi*lambda*dnplus1^2,1+alphan/2,'upper') )*gamma(1+alphan/2);
                    Ucal_term = gamma_difference/Gn_vec(end-aa+1)/(pi*lambda)^(alphan/2);
                    Ucal_sum = Ucal_sum + Ucal_term;
                end
                Ucal = P0/mu_UE*Ucal_sum;
                
            end
            %Effective transmit power per UE coefficient in W
            ETP_total =  K*Ucal*(tau_p + (tau_c-tau_p)*(ULfraction))/tau_c;
            
            %Total consumed power
            PC_total(m,k) = ETP_total + P_CP;
            
            %Compute EE with RZF
            EE_UL(m,k) = ThroughputPerCell(m,k)./PC_total(m,k);
            
        end
        
    end
    
end

%% Plot simulation results
[Krange_plot, Mrange_plot] = meshgrid(Krange,Mrange);
%Plot Figure 5.14a
figure;
hold on; box on; grid on;
surfc(Krange_plot, Mrange_plot, EE_UL/1e6, 'LineStyle', 'none')
colormap(autumn)
hold on
contour3(Krange_plot,Mrange_plot,EE_UL/1e6,10,'k')
[eeopt,I] = max(EE_UL(:));
[row,col] = ind2sub(size(EE_UL),I); % 2D maximizer
kopt = Krange_plot(1,col);  % Nms maximizer(Users=fixed)
mopt = Mrange_plot(row,1);  % Nbs maximizer (Users=fixed)
hold on
plot3(kopt,mopt,EE_UL(row,col)/1e6,'k^','MarkerSize',22,'MarkerFaceColor','black');
hold on
plot3(kopt,mopt,min(min(EE_UL)/1e6),'k^','MarkerSize',22,'MarkerFaceColor','black');

view([-37 30]);
xlabel('Number of UEs (K)')
ylabel('Number of BS antennas (M)')
zlabel('EE [Mbit/Joule]')
% xlabel('Number of UEs (K)')
% ylabel('Number of BS antennas (M)');
% zlabel('Energy efficiency [Mbit/Joule]');
% title('EE(M,K) with ZF combiner')
set(gca,'FontSize',20)
% Set Y-axis ticks every 50 units
yticks(0:50:250);

%% ALTERNATING OPTIMIZATION
if flagcombScheme == 2  %ZF combining
    %Initial point of the algorithm
    Mstar0 = 200;
    Kstar0 = 3;
    
    %Initialization
    iterMax = 10;
    Kstar_vec = [Kstar0, NaN*ones(1,iterMax)];
    Mstar_vec = [Mstar0, NaN*ones(1,iterMax)];
    EEstar_vec = [EE_UL(Mstar0,Kstar0), NaN*ones(1,iterMax)];
    zetastar_vec = [ZetaOpt(Mstar0,Kstar0), NaN*ones(1,iterMax)];
    
    %Plot the progress of the alternating optimization algorithm
    hold on; box on;
%     plot3(Kstar_vec(1),Mstar_vec(1),EEstar_vec(1)/1e6,'--b*','LineWidth',1,'MarkerFaceColor','blue','MarkerSize',20);
    %     plot3(Kstar_vec(2,1),Mstar_vec(2,1),EEstar_vec(2,1)/1e6,'--g*','LineWidth',3,'MarkerFaceColor','green','MarkerSize',20);
    
    %Compute CP calligraphic parameters
    [C0cal,C1cal,C2cal,C3cal,D0cal,D1cal,D2cal,E0cal] = functionCPparameters(tau_c,cp,Ucal);
    
    iter = 1;
    %Go through all iterations
    while iter<iterMax
        
        iter = iter + 1;
        
        %Extract the current values of K and M
        Kstar = Kstar_vec(iter-1);
        
        %% Step 3: Update M
        if pathloss == 1    %single-slope
            
            %Compute the parameters in Eqs. (31)-(36)
            a0 = gammaval*Kstar^2/tau_c*1/(alpha-1);
            a1 = (gammaval*Kstar/tau_c)*(Kstar*((2/(alpha-2))^2 + 2/(alpha-2) - 1/(alpha-1)) + 2/(alpha-2)/SNR);
            a2 = Kstar;
            a3 = Kstar*(1 + gammaval*(1/SNRp + 2/(alpha-2)*(1+1/SNRp))) + gammaval/SNR*(1+1/SNRp);
            a6 = C2cal*Kstar*tau_c;
            a5 = C0cal + C1cal*Kstar + C3cal*Kstar^3;
            a4 = D0cal*Kstar + D1cal*Kstar^2 + D2cal*Kstar^3;
            
        elseif pathloss == 2    %multislope
            
            %Compute the parameters in Eqs. (31)-(36)
            a0 = gammaval*Kstar^2/(tau_c*mu_kappa2);
            a1 = (gammaval*Kstar/tau_c)*(Kstar*(mu_kappa1^2 + mu_kappa1 - mu_kappa2) + mu_kappa1/SNR);
            a2 = Kstar;
            a3 = Kstar*(1 + gammaval*(1/SNRp + mu_kappa1*(1+1/SNRp))) + gammaval/SNR*(1+1/SNRp);
            a6 = C2cal*Kstar*tau_c;
            a5 = C0cal + C1cal*Kstar + C3cal*Kstar^3;
            a4 = D0cal*Kstar + D1cal*Kstar^2 + D2cal*Kstar^3;
            
        end
        
        r0 = a2-a0; r1 = a1+a3;
        q2 = a2*a4; q1 = -(a2*a5 - a3*a4 - a0*a6); q0 = (a1*a6 + a3*a5);
        %
        cbar1 = ( r1/r0 + sqrt(-q0/q2 - q1/q2*r1/r0 + (r1/r0)^2) );
        cbar2 = ( r1/r0 - sqrt(-q0/q2 - q1/q2*r1/r0 + (r1/r0)^2) );
        cbarStar_vec = max([cbar1,cbar2],[],2);
        % Checking feasibility
        cbar_c1 = (a1 + a3)/(a2 - a0);
        rho_c2 = Kstar/tau_c;
        cbar_c2 = -(a1 + rho_c2*a3)/(a0 - rho_c2*a2);
        %
        cbarStar = max([cbarStar_vec,cbar_c1,cbar_c2],[],2);
        
        %Compute the resulting M value
        Mstar = round(cbarStar*Kstar);
        
        %% Step 3: Update K
        if pathloss == 1    %single-slope
            
            %Compute the parameters in Eqs. (31)-(36)
            b0 = ((2/(alpha-2))^2 + cbarStar*1/(alpha-1) + 2/(alpha-2) - 1/(alpha-1))*gammaval/tau_c;
            b1 = 2/(alpha-2)*gammaval/SNR/tau_c;
            b2 = cbarStar - 1 - gammaval*2/(alpha-2)*(1+1/SNRp) - gammaval/SNRp;
            b3 = gammaval/SNR*(1+1/SNRp);
            b4 = C0cal;
            b5 = C1cal + D0cal*cbarStar;
            b6 = D1cal*cbarStar;
            b7 = C3cal + D2cal*cbarStar;
            b8 = tau_c*C2cal;
            
        elseif pathloss == 2    %multislope
            
            %Compute the parameters in Eqs. (31)-(36)
            b0 = (mu_kappa1^2 + cbarStar*mu_kappa2 + mu_kappa1 - mu_kappa2)*gammaval/tau_c;
            b1 = mu_kappa1*gammaval/SNR/tau_c;
            b2 = cbarStar - 1 - gammaval*mu_kappa1*(1+1/SNRp) - gammaval/SNRp;
            b3 = gammaval/SNR*(1+1/SNRp);
            b4 = C0cal;
            b5 = C1cal + D0cal*cbarStar;
            b6 = D1cal*cbarStar;
            b7 = C3cal + D2cal*cbarStar;
            b8 = tau_c*C2cal;
            
        end
        
        if flagOptK == 1    %exhaustive search
            
            m1 = b3;  m2 = (b2-b1);  m3 = b0;
            n0 = b3*b4;
            n1 = b2*b4 - b3*b5;
            n2 = b2*b5 - b3*b6 - b1*b8;
            n3 = b2*b6 - b3*b7 - b0*b8;
            n4 = b2*b7;
            
            Kcandidates = Krange;
            num = -m3*Kcandidates.^3 + m2*Kcandidates.^2 - m1*Kcandidates;
            den = n4*Kcandidates.^4 + n3*Kcandidates.^3 + n2*Kcandidates.^2 + n1*Kcandidates - n0;
            EE_Kcandidates = num./den;
            
            % Checking feasibility
            % Checking feasibility
            K_c1u = max([(-1+sqrt(1 - 4*b0*b3/(b1-b2)^2))*(b1-b2)/2/b0  (-1-sqrt(1 - 4*b0*b3/(b1-b2)^2))*(b1-b2)/2/b0]);
            K_c1d = min([(-1+sqrt(1 - 4*b0*b3/(b1-b2)^2))*(b1-b2)/2/b0  (-1-sqrt(1 - 4*b0*b3/(b1-b2)^2))*(b1-b2)/2/b0]);
            K_c2 = (b1*tau_c + b3)/(b2 - b0*tau_c);
            K_c3 = Mstar_vec(1,1);
            %
            EE_Kcandidates(Kcandidates<K_c2) = 0;
            EE_Kcandidates(Kcandidates<K_c1d | Kcandidates>K_c1u) = 0;
            EE_Kcandidates(Kcandidates>K_c3) = 0;
            %
            [EE_Kstar,Kstar] = max(EE_Kcandidates);
            
        elseif flagOptK == 2    %simplified
            
            Kstar = round(C0cal/(C1cal + D0cal*cbarStar)*(sqrt(1+(b2-b1)/b0*(C1cal + D0cal*cbarStar)/C0cal) - 1)) ;
            
        end
        
        %% Joint update M and K
        %Extract the current K and M values, which are used to compute the current EE value
        Kstar_vec(iter) = Kstar;
        Mstar_vec(iter) = Mstar;
        
        %% Compute optimum point
        %Compute the current EE value
        EEstar_vec(iter) = EE_UL(Mstar_vec(iter),Kstar_vec(iter));
        %Compute the current pilot reuse factor
        zetastar_vec(iter) = ZetaOpt(Mstar,Kstar);
        
        %Plot the progress of the alternating optimization algorithm
        hold on; box on;
        plot3([Kstar_vec(iter-1) Kstar_vec(iter)],[Mstar_vec(iter-1) Mstar_vec(iter)],[EEstar_vec(iter-1) EEstar_vec(iter)]/1e6,'-bp','LineWidth',3,'MarkerFaceColor','blue','MarkerSize',24);
        
        %% Exit condition
        %Evaluate if the algorithm has reached a steady point
        cond1 = Mstar_vec(iter)==Mstar_vec(iter-1) && Kstar_vec(iter)==Kstar_vec(iter-1);
        %         cond2 = Mstar_vec(2,iter)==Mstar_vec(2,iter-1) && Kstar_vec(2,iter)==Kstar_vec(2,iter-1);
        if cond1
            break;
        end
        
        if iter==iterMax
            display('Maximum number of iteration reached')
        end
        
    end
    
end

%Display the optimal network deployment setup
SINR = gammaval
AvgRate = rateval
EEopt = eeopt/1e6
ASEopt = lambda*ThroughputPerCell(mopt,kopt)/1e6
APCopt = lambda*PC_total(mopt,kopt)
Mopt = mopt
Kopt = kopt
Zetaopt = ZetaOpt(mopt,kopt)
Reuseopt = 1/ZetaOpt(mopt,kopt)*100


