function [P_MR,P_ZF,P_MMMSE,cp,P_RZF,P_SMMSE] = functionCPcomputationvsLvsP(lambdaRange,fRange,M,K,L,B,tau_c,tau_p,ULfraction,valueset,sumSE_MR_UL,sumSE_ZF_UL,sumSE_MMMSE_UL,sumSE_RZF_UL,sumSE_SMMSE_UL)
 
%This Matlab function is called by the script plot_EEvsThroughput.m in the article:
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
%
%This function returns the total CP with different processing schemes,
%using the model defined Section 2, the parameters in Table 1 and the results 
% in Table 2
%
%
%INPUT:
%lambdaRange    = Vector with different values of BS densities
%fRange         = Vector with different values of pilot reuse factors
%M              = Number of BS antennas
%K              = Number of UEs per cell
%L              = Number of cells
%sumSE_MR_UL    = Average sum uplink SE per cell with MR (same size as Mrange)
%sumSE_RZF_UL   = Average sum uplink SE per cell with RZF (same size as Mrange)
%sumSE_MMMSE_UL = Average sum uplink SE per cell with M-MMSE (same size as Mrange)
%sumSE_ZF_UL    = Average sum uplink SE per cell with ZF (same size as Mrange)
%sumSE_SMMSE_UL = Average sum uplink SE per cell with S-MMSE (same size as Mrange)
%B              = Bandwidth
%tau_c          = Length of coherence block
%tau_p          = Length of pilot sequence
%
%valueset       = Select which one of the value sets in Table 5.3 that is
%               considered. Either valueset=1 or valueset=2 are possible.
%ULfraction     = Fractions of UL and DL data - (tau_c - tau_p)/(tau_c) accounts for the
%               remaining samples after pilot transmission.
%OUTPUT:
%P_MR    = Vector of same size as Mrange with total CP for MR processing
%P_RZF   = Same as P_MR but with RZF processing
%P_MMMSE = Same as P_MMMSE but with M-MMSE processing
%P_ZF    = Same as P_ZF but with XF processing
%P_SMMSE = Same as P_SMMSE but with S-MMSE processing

%Obtain CP model coefficients for one of the value sets in Table 5.3
[P_FIX,P_LO,P_BS,P_UE,P_COD,P_DEC,L_BS,P_BT] = functionCPmodel(valueset);

%Prepare to store simulation results
P_MR = zeros(length(lambdaRange),length(fRange));
P_TC = zeros(length(lambdaRange),length(fRange));
P_CE = zeros(length(lambdaRange),length(fRange));
P_SP_RT = zeros(length(lambdaRange),length(fRange));
% P_SP_DL = zeros(length(lambdaRange),length(fRange));
P_SP_UL_MR = zeros(length(lambdaRange),length(fRange));
P_SAME = zeros(length(lambdaRange),length(fRange));

%Compute CP for coding and decoding using (5.35)
P_CD_MR = (P_COD + P_DEC)*B*sumSE_MR_UL;

%Compute CP for backhaul traffic using (5.36)
P_BH_MR = P_BT*B*sumSE_MR_UL;

%Repeat computations for ZF
if nargin>11
    P_ZF = zeros(length(lambdaRange),length(fRange));
    P_CD_ZF = (P_COD + P_DEC)*B*sumSE_ZF_UL;
    P_BH_ZF = P_BT*B*sumSE_ZF_UL;
    P_SP_UL_ZF = zeros(length(lambdaRange),length(fRange));
end

%Repeat computations for M-MMSE
if nargin>12
    P_MMMSE = zeros(length(lambdaRange),length(fRange));
    P_CD_MMMSE = (P_COD + P_DEC)*B*sumSE_MMMSE_UL;
    P_BH_MMMSE = P_BT*B*sumSE_MMMSE_UL;
    P_SP_UL_MMMSE = zeros(length(lambdaRange),length(fRange));
end

%Repeat computations for RZF
if nargin>13
    P_RZF = zeros(length(lambdaRange),length(fRange));
    P_CD_RZF = (P_COD + P_DEC)*B*sumSE_RZF_UL;
    P_BH_RZF = P_BT*B*sumSE_RZF_UL;
    P_SP_UL_RZF = zeros(length(lambdaRange),length(fRange));
end

%Repeat computations for S-MMSE
if nargin>14
    P_SMMSE = zeros(length(lambdaRange),length(fRange));
    P_CD_SMMSE = (P_COD + P_DEC)*B*sumSE_SMMSE_UL;
    P_BH_SMMSE = P_BT*B*sumSE_SMMSE_UL;
    P_SP_UL_SMMSE = zeros(length(lambdaRange),length(fRange));
end


% Go through all number of BS densities
for l = 1:length(lambdaRange)
    
    %Go through all number of Pilot reuse factors
    for f = 1:length(fRange)
        
        %Compute CP for transceiver chains using (5.34)
        P_TC(l,f) = M*P_BS + P_LO + K*P_UE;
        
        %Compute CP for channel estimation with all other schemes, where only
        %the channels to UEs in other cells are estimated, using (5.37)
        P_CE(l,f) = 3*K*B/(tau_c*L_BS)*(M*tau_p + M);  %3*K*B/(tau_c*L_BS)*(M*tau_p + M^2);
        
        %Compute CP for UL reception
        P_SP_RT(l,f) = ULfraction*(tau_c - tau_p)*3*B/(tau_c*L_BS)*M*K;
        
        %         %Compute CP for computation of precoding vectors
        %         P_SP_DL(l,f) = 4*M*K*B/(tau_c*L_BS);
        
        %Sum up the power terms that are independent of the processing scheme
        P_SAME(l,f) = P_FIX + P_TC(l,f) + P_SP_RT(l,f); %+ P_SP_DL(l,f);
        
        
        %Compute CP for computation of the combining vectors with different
        %schemes, based on Table 5.2
        P_SP_UL_MR(l,f) = 7*B*K/(tau_c*L_BS);
        
        %Compute the final CP values
        P_MR(l,f) = P_SAME(l,f) + P_CE(l,f) + P_CD_MR(l,f) + P_BH_MR(l,f) + P_SP_UL_MR(l,f);
        
        %Repeat same computations for ZF
        if nargin>11
            P_SP_UL_ZF(l,f) = 3*B*(3*K^2*M/2 + M*K/2 + (K^3 - K)/3 + (7/3)*K)/(tau_c*L_BS);
            P_ZF(l,f) = P_SAME(l,f) + P_CE(l,f) + P_CD_ZF(l,f) + P_BH_ZF(l,f) + P_SP_UL_ZF(l,f);
        end
        
        
        %Repeat same computations for M-MMSE
        if nargin>12
            P_SP_UL_MMMSE(l,f) = 3*B*(L*(M^2 + 3*M)*K/2 + M^3/3 + (M^2 - M)*K + 2*M + M*tau_p*(tau_p-K))/(tau_c*L_BS);  % 3*B*(L*(3*M^2 + M)*K/2 + M^3/3 + 2*M + M*tau_p*(tau_p-K))/(tau_c*L_BS);
            P_MMMSE(l,f) = P_SAME(l,f) + P_CE(l,f) + P_CD_MMMSE(l,f) + P_BH_MMMSE(l,f) + P_SP_UL_MMMSE(l,f);
        end
        
        %Repeat same computations for RZF
        if nargin>13
            P_SP_UL_RZF(l,f) = 3*B*(3*K^2*M/2 + 3*M*K/2 + (K^3 - K)/3 + (7/3)*K)/(tau_c*L_BS);
            P_RZF(l,f) = P_SAME(l,f) + P_CE(l,f) + P_CD_RZF(l,f) + P_BH_RZF(l,f) + P_SP_UL_RZF(l,f);
        end
        
        %Repeat same computations for S-MMSE
        if nargin>14
            P_SP_UL_SMMSE(l,f) = 3*B*(3*M^2*K/2 + M*K/2 + (M^3 - M)/3 + (7/3)*M)/(tau_c*L_BS);
            P_SMMSE(l,f) = P_SAME(l,f) + P_CE(l,f) + P_CD_SMMSE(l,f) + P_BH_SMMSE(l,f) + P_SP_UL_SMMSE(l,f);
        end
        
    end
    
end

%Store the Circuit Power consumption terms
cp.P_SAME = P_SAME;
cp.P_CE = P_CE;
cp.P_CD_MR = P_CD_MR;
cp.P_BH_MR = P_BH_MR;
cp.P_SP_UL_MR = P_SP_UL_MR;

end
