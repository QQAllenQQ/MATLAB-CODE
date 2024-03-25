function [C0cal,C1cal,C2cal,C3cal,D0cal,D1cal,D2cal,E0cal] = functionCPparameters(tau_c,cp,ucal)

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
%This function returns the Auxiliary power consumption coefficient used in
%Section 5A.
%
%INPUT:
%tau_c      = Length of coherence block
%cp         = Matlab structure including Circuit Power consumption terms
%ucal       = Auxiliary power coefficient given in Corollary 3
%
%OUTPUT:
%C0cal      = PC coefficient associated to constant term 
%C1cal      = PC coefficient scaling with K
%C2cal      = PC coefficient scaling with K^2
%C3cal      = PC coefficient scaling with K^3
%D0cal      = PC coefficient scaling with M
%D1cal      = PC coefficient scaling with (K x M)
%D2cal      = PC coefficient scaling with (K^2 x M)
%E0cal      = PC coefficient multiplying the Area Rate

%Load CP model
L_BS = cp.L_BS;
P_FIX = cp.P_FIX;
P_LO = cp.P_LO;
P_BS = cp.P_BS;
P_UE = cp.P_UE;
P_COD = cp.P_COD;
P_DEC = cp.P_DEC;
P_BT = cp.P_BT;

%Compute CP calligraphic parameters
aux = 3/tau_c/L_BS;
C0cal = (P_FIX + P_LO);
C1calbar = (P_UE + aux*(2-1/3));
C3cal = (aux/3);
D0cal = (P_BS);
D1cal = (aux*(1/2 + tau_c + 2));
D2cal = (aux*3/2);
E0cal = (P_COD + P_DEC + P_BT);

Ctx1 = (ucal*(1+1/tau_c));
C1cal = C1calbar + Ctx1;
C2cal = (ucal/tau_c);


end