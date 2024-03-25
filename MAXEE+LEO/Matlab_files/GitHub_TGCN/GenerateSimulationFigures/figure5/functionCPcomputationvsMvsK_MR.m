function [P_MR,cp] = functionCPcomputationvsMvsK_MR(M,K,B,tau_c,tau_p,ULfraction,valueset,ThroughputPerCell)

%This Matlab function is called by the script plot_EEvsMvsK.m in the article:
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
%This function returns the total CP for MR processing scheme only,
%using the model defined Section 2, the parameters in Table 1 and the results 
% in Table 2
%
%
%INPUT:
%M                  = Number of BS antennas
%K                  = Number of UEs per cell
%B                  = Bandwidth
%tau_c              = Length of coherence block
%tau_p              = Length of pilot sequence
%ULfraction         = Fractions of UL and DL data - (tau_c - tau_p)/(tau_c) accounts for the
%                   remaining samples after pilot transmission
%valueset           = Select the value set of the power coefficients that is considered 
%                   (set valueset=2 to use the hardware parameters listen in Table 1 and 
%                   generate the results as in the article)
%ThroughputPerCell  = Uplink average ergodic Area Rate 
%
%OUTPUT:
%P_MR    = Vector of same size as Mrange with total CP for MR processing
%

%Obtain CP model coefficients for one of the value sets in Table 5.3
[P_FIX,P_LO,P_BS,P_UE,P_COD,P_DEC,L_BS,P_BT] = functionCPmodel(valueset);

%Coding/Decoding and Backhauling costs
P_CD_MR = (P_COD + P_DEC)*ThroughputPerCell;
P_BH_MR = P_BT*ThroughputPerCell;

%Compute CP for transceiver chains using (5.34)
P_TC = M*P_BS + P_LO + K*P_UE;

%Compute CP for channel estimation with all other schemes, where only
%the channels to UEs in other cells are estimated, using (5.37)
P_CE = 3*K*B/(tau_c*L_BS)*(M*tau_p + M);  %3*K*B/(tau_c*L_BS)*(M*tau_p + M^2);

%Compute CP for UL reception
P_SP_RT = ULfraction*(tau_c - tau_p)*3*B/(tau_c*L_BS)*M*K;

%         %Compute CP for computation of precoding vectors
%         P_SP_DL(m,k) = 4*M*K*B/(tau_c*L_BS);

%Sum up the power terms that are independent of the processing scheme
P_SAME = P_FIX + P_TC + P_SP_RT; %+ P_SP_DL(m,k);

%Compute CP for computation of the combining vectors with different
%schemes, based on Table 5.2
P_SP_UL_MR = 7*B*K/(tau_c*L_BS);

%Compute the final CP values        
P_MR = P_SAME + P_CE + P_CD_MR + P_BH_MR + P_SP_UL_MR;

%Store the Circuit Power consumption terms
cp.L_BS = L_BS;
cp.P_FIX = P_FIX;
cp.P_LO = P_LO;
cp.P_BS = P_BS;
cp.P_UE = P_UE;
cp.P_COD = P_COD;
cp.P_DEC = P_DEC;
cp.P_BT = P_BT;

end
