%This Matlab script is called by the script plot_EEvsMvsK.m in the article:
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
%This script provides the script plot_EEvsMvsK.m with the nedeed paramters.
%

%Select the value set of the power coefficients that is considered 
%(set valueset=2 to use the hardware parameters listen in Table 1 and 
% generate the results as in the article)
valueset = 2;

%Select the BS densities
lambda = 10;

% Number of BSs
L = lambda*(squareLength/1000)^2;

%Select the values of log2(1+gamma) that should be considered (\gamma > \alpha-1, for feasibility)
gammaval = 7;
rateval = log2(1 + gammaval);

%Pathloss exponent
alpha = 3.76;

%path-loss exponents used in the multi slope model
pathLossExp_vec = Pathloss_model.pathLossExp_vec;

%reference distances at which a change in power decandence occur used in the multi slope model
distancesPathLossExp_vec = Pathloss_model.distancesPathLossExp_vec;

%multislope pathloss coefficients
Gn_vec = Pathloss_model.Gn_vec;

% %Number of different path-loss exponents used in the multi slope model
% nbrOfSlopes = 3;
% %path-loss exponents used in the multi slope model
% pathLossExp_vec = linspace(3.76,2.1,nbrOfSlopes);
% distancesPathLossExp_vec = linspace(0,1000,nbrOfSlopes+1);
% %Multi slope (Andrea)
% Gn_vec = [Pathlossat1km, distancesPathLossExp_vec(end).^-pathLossExp_vec(2:nbrOfSlopes)];

