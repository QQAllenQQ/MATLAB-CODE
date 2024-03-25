%This Matlab script is called by the script main.m in the article:
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
%This script includes the channel propagation and power allocation
%parameters that are used in the main code.


%Set the length in meters of the total square area coverage
squareLength = 1000; % Area coverage = (squareLength/1000)^2 [Km^2]

%Communication bandwidth
B = 20e6;

%Select length of coherence block
tau_c = 200;

%Fractions of UL and DL data - (tau_c - tau_p)/(tau_c) accounts for the
%remaining samples after pilot transmission
ULfraction = 1/3;

%Use the approximation of the Gaussian local scattering model
accuracy = 2;

%Angular standard deviation in the local scattering model (in degrees)
ASDdeg = 10;

%Total uplink transmit power per UE (mW)
p = 100;

%Noise figure at the BS (in dB)
noiseFigure = 7;

%Compute noise power
noiseVariancedBm = -174 + 10*log10(B) + noiseFigure;

%Set the SNR0 in dB in the inverse statistical power control policy
%(payload data)
DeltadB = 5;

%Set the SNRp in dB in the inverse statistical power control policy
%(pilot signaling)
DeltadB_pilot = 15;

%Total uplink transmit power per UE (W)
P0dBm = DeltadB + noiseVariancedBm;

%Total uplink transmit power per UE (W)
P0 = 10^(0.1*(P0dBm - 30));

%PA efficiency UEs and BSs
mu_UE = 0.4;
mu_BS = 0.5;

%% Pathloss model
%Standard deviation of shadow fading
Pathloss_model.sigma_sf = 0;

%Minimum distance between BSs and UEs
Pathloss_model.minDistance = 35;

%Define the antenna spacing (in number of wavelengths)
Pathloss_model.antennaSpacing = 1/2; %Half wavelength distance

%%% Single-slope
%Pathloss exponent
alpha = 3.76;
Pathloss_model.alpha = alpha;

%Pathloss at reference distance of 1 km
freq_MHz = 2e3; %carrier frequency in MHz
height_BS = 10; %BS antenna height in m
height_UE = 1.65;   %UE antenna height in m
Pathloss_model.Pathlossat1km_dB = 46.3 + 33.9*log10(freq_MHz) - 13.82*log10(height_BS) - (1.1*log10(freq_MHz)-0.7)*height_UE + (1.56*log10(freq_MHz)-0.8);
% Pathlossat1km_dB = 148.1;
Pathlossat1km = 10^(-Pathloss_model.Pathlossat1km_dB/10);

%Average channel gain in dB at a reference distance of 1 meter. Note that
%-35.3 dB corresponds to -148.1 dB at 1 km, using pathloss exponent 3.76
Pathloss_model.constantTerm = -31.1;

%%% Multislope
%Number of different path-loss exponents used in the multi slope model
Pathloss_model.nbrOfSlopes = 3;
% Pathloss_model.nbrOfSlopes = 4;

%path-loss exponents used in the multi slope model
Pathloss_model.pathLossExp_vec = [4, 2.01, 0];
% exp0 = 0;   exp1 = 2.01;   exp2 = 3;   exp3 = 4;
% Pathloss_model.pathLossExp_vec = [exp3 exp2 exp1 exp0];

%reference distances at which a change in power decandence occur used in the multi slope model
Rcutoff = 4*height_BS*height_UE/3e2*freq_MHz;
Pathloss_model.distancesPathLossExp_vec = [0, 10, Rcutoff, 1000];
% d0 = 1; d1 = 20;    d2 = 200;
% Pathloss_model.distancesPathLossExp_vec = [0 d0 d1 d2 1000];

%multislope pathloss coefficients
Pathloss_model.Gn_vec = [Pathlossat1km, Pathloss_model.distancesPathLossExp_vec(end).^-Pathloss_model.pathLossExp_vec(2:Pathloss_model.nbrOfSlopes)];
% dref = 1e3;
% startpoint = 10^(-4);
% const = d1^(exp2)*d2^(exp3);
% Pathloss_model.Gn_vec = startpoint*[dref^(-exp3)*const*d1^(-exp1)*d2^(-exp2), dref^(-exp2)*const*d1^(-exp1)*d2^(-exp3), ...
%     dref^(-exp1)*const*d1^(-exp2)*d2^(-exp3), 1*dref^(-exp0)];

