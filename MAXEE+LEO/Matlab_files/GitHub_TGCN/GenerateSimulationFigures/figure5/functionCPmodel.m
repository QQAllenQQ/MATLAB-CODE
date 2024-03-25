function [P_FIX,P_LO,P_BS,P_UE,P_COD,P_DEC,L_BS,P_BT] = functionCPmodel(valueset)

%This Matlab function is called by the function functionCPcomputationvsMvsK_ZF.m in the article:
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
%This function returns the values of the CP model coefficients in Table 1.
%
%
%INPUT:
%valueset = Select which one of the value sets in Table 5.3 that is
%           considered. Either valueset=1 or valueset=2 are possible.
%
%OUTPUT:
% P_FIX = Fixed circuit power per BS (in Watt)
% P_LO  = Circuit power per LO (in Watt)
% P_BS  = Circuit power per BS antenna (in Watt)
% P_UE  = Circuit power UE (in Watt)
% P_COD = Circuit power for encoding
% P_DEC = Circuit power for decoding
% L_BS  = BS computational efficiency
% P_BT  = Circuit power for backhaul


%% Define parameter values for Value set 1
if valueset == 1
    
    P_FIX = 10;
    
    P_LO = 0.2;
    
    P_BS = 0.4;
    
    P_UE = 0.2;
    
    P_COD = 0.1*10^(-9);
    
    P_DEC = 0.8*10^(-9);
    
    L_BS = 75*10^9;
    
    P_BT = 0.25*10^(-9);
    
    % Define parameter values for Value set 2
elseif valueset == 2
    
    P_FIX = 5;
    
    P_LO = 0.1;
    
    P_BS = 0.2;
    
    P_UE = 0.1;
    
    P_COD = 0.01*10^(-9);
    
    P_DEC = 0.08*10^(-9);
    
    L_BS = 750*10^9;
    
    P_BT = 0.025*10^(-9);
    
end


end
