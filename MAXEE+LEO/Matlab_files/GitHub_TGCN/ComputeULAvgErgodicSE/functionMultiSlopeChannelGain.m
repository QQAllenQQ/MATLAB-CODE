function [channelGaindB] = functionMultiSlopeChannelGain(distancesBSj,pathLossExp_vec,distancesPathLossExp_vec,Gn_vec)

%This Matlab function is called by the function functionExampleSetup.m in the article:
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
%This functions computes the average channel gain of the channel between UEs at random
%locations and the BSs as defined in Section 2.
%
%
%INPUT:
%distancesBSj               = K x 1 vector including the distances from all the 
%                           UEs in cell l to BS j with a wrap around topology, 
%                           where the shortest distance between a UE and the 
%                           nine different locations of a BS is considered
%pathLossExp_vec            = Path loss exponents used in the multi slope model
%distancesPathLossExp_vec   = Reference distances at which a change in power 
%                           decandence occur used in the multi slope model
%Gn_vec                     = multislope pathloss coefficients 
%
%OUTPUT:
%channelGaindB              = K x 1 vector including the average channel
%                           gain of all the UE in cell l to BS j


%Retrieve the number of UEs K and the number of slopes in the multi-slope
%model
K = length(distancesBSj);
nbrOfSlopes = length(pathLossExp_vec);

% %Compute the large scale fading coefficients (assign path-loss exponent according to the distances BS_l - UE_l,j)
% channelGain = NaN*ones(K,1);
% indexesPathLoss = distancesBSj<distancesPathLossExp_vec(2);
% channelGain(indexesPathLoss) = Kn_vec(1)*(distancesBSj(indexesPathLoss).^-pathLossExp_vec(1));
% for ii = 2:nbrOfSlopes-1
%     indexesPathLossii = distancesBSj>=distancesPathLossExp_vec(ii)&distancesBSj<distancesPathLossExp_vec(ii+1);
%     channelGain(indexesPathLossii) = Kn_vec(ii)*(distancesBSj(indexesPathLossii).^-pathLossExp_vec(ii));
% end
% channelGain(distancesBSj>=distancesPathLossExp_vec(end-1)) = Kn_vec(end)*(distancesBSj(distancesBSj>=distancesPathLossExp_vec(end-1)).^-pathLossExp_vec(end-1));

%Compute the large-scale fading coefficients (assign path-loss exponent according to the distances BS_l - UE_l,j)
channelGain = NaN*ones(K,1);
for ii = 1:nbrOfSlopes
    indexesPathLossii = find(distancesBSj>=distancesPathLossExp_vec(ii)&distancesBSj<distancesPathLossExp_vec(ii+1));
    channelGain(indexesPathLossii) = Gn_vec(end-ii+1)*((distancesBSj(indexesPathLossii)/1e3).^-pathLossExp_vec(end-ii+1));
end

%Conversion to db
channelGaindB = 10*log10(channelGain);

end

