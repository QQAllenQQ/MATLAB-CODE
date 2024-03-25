function [R,channelGaindB] = functionExampleSetup(L,K,M,accuracy,ASDdeg,squareLength,Pathloss_model,deployment,pathloss)

%This Matlab function is called by the script main.m in the article:
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
%This function generates the channel statistics (large scale fading) between UEs at random
%locations and the BSs as defined in Section 2.
%
%Note that all distances in this code are measured in meters.
%
%INPUT:
%L                  = Number of BSs and cells
%K                  = Number of UEs per cell
%M                  = Number of antennas per BS
%accuracy           = Compute exact correlation matrices from the local scattering
%                   model if approx=1. Compute a small-angle approximation of the
%                   model if approx=2
%ASDdeg             = Angular standard deviation around the nominal angle
%                   (measured in degrees)
%squareLength       = Length in meters of the total square area coverage
%Pathloss_model     = Matlab structure containing all the path loss
%                   parameters for the multi slope model (e.g., number of slopes, path loss exponents)
%deployment         = Generates BS locations over (1) a fixed grid or (2) 
%                   as a Homogeneous Poisson Point Process H-PPP
%pathloss           = Generates the large scale path-loss coefficient using (1) single slope or
%                   (2) multi slope model
%
%OUTPUT:
%R             = M x M x K x L x L matrix with spatial correlation matrices
%                for all UEs in the network. R(:,:,k,j,l) is the correlation
%                matrix for the channel between UE k in cell j and the BS
%                in cell l. This matrix is normalized such that trace(R)=M.
%channelGaindB = K x L x L matrix containing the average channel gains in
%                dB of all the channels. The product
%                R(:,:,k,j,l)*10^(channelGaindB(k,j,l)/10) is the full
%                spatial channel correlation matrix.


%% Model parameters

%Standard deviation of shadow fading
sigma_sf = Pathloss_model.sigma_sf;

%Minimum distance between BSs and UEs
minDistance = Pathloss_model.minDistance;

%Define the antenna spacing (in number of wavelengths)
antennaSpacing = Pathloss_model.antennaSpacing; %Half wavelength distance

if pathloss == 1
    
    %Pathloss exponent
    alpha = Pathloss_model.alpha;
    
    %Average channel gain in dB at a reference distance of 1 meter. Note that
    %-35.3 dB corresponds to -148.1 dB at 1 km, using pathloss exponent 3.76
    constantTerm = Pathloss_model.constantTerm;
    
elseif pathloss == 2
    
%     %Number of different path-loss exponents used in the multi slope model
%     nbrOfSlopes = 6;
%     %Multi-slope model Pathloss exponent
%     pathLossExp_vec = linspace(2,3.76,nbrOfSlopes);
%     distancesPathLossExp_vec = linspace(0,1000,nbrOfSlopes+1);
%     %Multi slope Andrews
%     DeltaPathLossExp = pathLossExp_vec(2)-pathLossExp_vec(1);
%     Risum = log(prod(distancesPathLossExp_vec(2:end-2)));
%     omegaPathLossExp = DeltaPathLossExp*Risum / log(1/Pathlossat1km*(distancesPathLossExp_vec(end))^-pathLossExp_vec(end-1));
%     Kn_vec = zeros(1,nbrOfSlopes);
%     Kn_vec(1) = 1;
%     for indexa = 2:nbrOfSlopes-1
%         Kn_vec(indexa) = Kn_vec(indexa-1)*(distancesPathLossExp_vec(indexa)).^-((pathLossExp_vec(indexa)-pathLossExp_vec(indexa-1))/omegaPathLossExp);
%     end
%     Kn_vec(end) = Kn_vec(end-1);

    %path-loss exponents used in the multi slope model
    pathLossExp_vec = Pathloss_model.pathLossExp_vec;
    
    %reference distances at which a change in power decandence occur used in the multi slope model
    distancesPathLossExp_vec = Pathloss_model.distancesPathLossExp_vec;
    
    %multislope pathloss coefficients
    Gn_vec = Pathloss_model.Gn_vec;
    
end

if deployment == 1
    
    %Number of BSs per dimension
    nbrBSsPerDim = sqrt(L);
    
    %Distance between BSs in vertical/horizontal direction
    interBSDistance = squareLength/nbrBSsPerDim;
    
    %Deploy BSs on the grid
    locationsGridHorizontal = repmat(interBSDistance/2:interBSDistance:squareLength-interBSDistance/2,[nbrBSsPerDim 1]);
    locationsGridVertical = locationsGridHorizontal';
    BSpositions = locationsGridHorizontal(:) + 1i*locationsGridVertical(:);
    
    %Compute all nine alternatives of the BS locations when using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);
    
elseif  deployment == 2
    
    %Random BS locations with uniform distribution
    BSpositions = (rand(L,1) + 1i*rand(L,1)) * squareLength;
    
    %Compute all nine alternatives of the BS locations when using wrap around
    wrapHorizontal = repmat([-squareLength 0 squareLength],[3 1]);
    wrapVertical = wrapHorizontal';
    wrapLocations = wrapHorizontal(:)' + 1i*wrapVertical(:)';
    BSpositionsWrapped = repmat(BSpositions,[1 length(wrapLocations)]) + repmat(wrapLocations,[L 1]);
    
    %Compute a size of the distance from a BS where its users need to be
    % contained.
    distances = abs(BSpositions);
    sortedBSdistances = sort(distances,'ascend');
    maxDistance = 2*sortedBSdistances(end);
    
end

%Prepare to put out UEs in the cells
UEpositions = zeros(K,L);
perBS = zeros(L,1);

%Prepare to store normalized spatial correlation matrices
R = zeros(M,M,K,L,L,length(ASDdeg));

%Prepare to store average channel gain numbers (in dB)
channelGaindB = zeros(K,L,L);


%% Go through all the cells
for l = 1:L
    
    if deployment == 1
        %Put out K UEs in the cell, uniformly at random. The procedure is
        %iterative since UEs that do not satisfy the minimum distance are
        %replaced with new UEs
        while perBS(l)<K
            
            %Put out new UEs
            UEremaining = K-perBS(l);
            posX = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
            posY = rand(UEremaining,1)*interBSDistance - interBSDistance/2;
            posXY = posX + 1i*posY;
            
            %Keep those that satisfy the minimum distance
            posXY = posXY(abs(posXY)>=minDistance);
            
            %Store new UEs
            UEpositions(perBS(l)+1:perBS(l)+length(posXY),l) = posXY + BSpositions(l);
            perBS(l) = perBS(l)+length(posXY);
            
        end
        
    elseif  deployment == 2
        
        %Put out K UEs in the cell, uniformly at random. The procedure is
        %iterative since UEs that do not satisfy the minimum distance are
        %replaced with new UEs
        while perBS(l)<K
            
            %Put out new UEs
            UEremaining = K-perBS(l);
            posX = rand(UEremaining,1)*maxDistance + real(BSpositions(l)) - maxDistance/2;
            posY = rand(UEremaining,1)*maxDistance + imag(BSpositions(l)) - maxDistance/2;
            posXY = mod(posX,squareLength) + 1i*mod(posY,squareLength);
            
            %%UPDATED FROM HERE
            
            %Find closest BS (with wrap around)
            for k = 1:UEremaining
                
                distancesUEtoBSs = abs(BSpositionsWrapped - repmat(posXY(k),size(BSpositionsWrapped)));
                
                [~,index] = min(min(distancesUEtoBSs,[],2),[],1);
                
                if index == l
                    
                    UEpositions(perBS(l)+1,l) = posXY(k);
                    perBS(l) = perBS(l)+1;
                end
                
            end
            
            %%TO HERE
            
        end
        
    end
    
    
    %Go through all BSs
    for j = 1:L
        
        %Compute the distance from the UEs in cell l to BS j with a wrap
        %around topology, where the shortest distance between a UE and the
        %nine different locations of a BS is considered
        [distancesBSj,whichpos] = min(abs( repmat(UEpositions(:,l),[1 size(BSpositionsWrapped,2)]) - repmat(BSpositionsWrapped(j,:),[K 1]) ),[],2);
        
        %Compute average channel gain using the large-scale fading model in
        %(2.3), while neglecting the shadow fading
        if pathloss == 1
            
            channelGaindB(:,l,j) = constantTerm - alpha*10*log10(distancesBSj);
            
        elseif pathloss == 2
            
            channelGaindB(:,l,j) = functionMultiSlopeChannelGain(distancesBSj,pathLossExp_vec,distancesPathLossExp_vec,Gn_vec);
            
        end
        
        %Compute nominal angles between UE k in cell l and BS j, and
        %generate spatial correlation matrices for the channels using the
        %local scattering model
        for k = 1:K
            
            angleBSj = angle(UEpositions(k,l)-BSpositionsWrapped(j,whichpos(k)));
            
            if accuracy == 1 %Use the exact implementation of the local scattering model
                
                for spr = 1:length(ASDdeg)
                    
                    R(:,:,k,l,j,spr) = functionRlocalscattering(M,angleBSj,ASDdeg(spr),antennaSpacing);
                    
                end
                
            elseif accuracy == 2 %Use the approximate implementation of the local scattering model
                
                for spr = 1:length(ASDdeg)
                    
                    R(:,:,k,l,j,spr) = functionRlocalscatteringApprox(M,angleBSj,ASDdeg(spr),antennaSpacing);
                    
                end
                
            end
            
        end
        
    end
    
    
    %Go through all UEs in cell l and generate shadow fading realizations
    for k = 1:K
        
        %Generate shadow fading realizations
        shadowing = sigma_sf*randn(1,1,L);
        channelGainShadowing = channelGaindB(k,l,:) + shadowing;
        
        %Check if another BS has a larger average channel gain to the UE
        %than BS l
        while channelGainShadowing(l) < max(channelGainShadowing)
            
            %Generate new shadow fading realizations (until all UEs in cell
            %l has its largest average channel gain from BS l)
            shadowing = sigma_sf*randn(1,1,L);
            channelGainShadowing = channelGaindB(k,l,:) + shadowing;
            
        end
        
        %Store average channel gains with shadowing fading
        channelGaindB(k,l,:) = channelGainShadowing;
        
    end
    
end

% %%
% %Generate random pilot pattern
% f = 4; % Pilot reuse factor
% pilotPattern = randi(f,L,1);
% % Cell array of colors
% % Colarray = rand(1,3*f);
% Colarray = [0.8980, 0.1686, 0.3137, 0.5882, 0.2941, 0, 0, 0, .5, 1, 0.6, 0];
% % Cell array of markers
% Markarray = ['p','*','o','d','^','x','>','<'];
% 
% figure;
% hold on; box on;
% voronoi(real(BSpositionsWrapped),imag(BSpositionsWrapped));
% plot(real(BSpositions),imag(BSpositions),'s','MarkerSize',1,'MarkerFaceColor','r');
% axis([0 squareLength 0 squareLength]);
% hold on;
% plot(real(UEpositions(1,1)),imag(UEpositions(1,1)),'^','MarkerSize',20,'color',Colarray(3*(pilotPattern(1)-1)+1:3*pilotPattern(1)),'Marker',Markarray(pilotPattern(1)),'MarkerFaceColor',Colarray(3*(pilotPattern(1)-1)+1:3*pilotPattern(1))');
% for l = 1:L
% %     pilotPattern(l)
% %     Carray(3*(pilotPattern(l)-1)+1:3*pilotPattern(l))
%     plot(real(UEpositions(:,l)),imag(UEpositions(:,l)),'*','MarkerSize',10,'color',Colarray(3*(pilotPattern(l)-1)+1:3*pilotPattern(l)),'Marker',Markarray(pilotPattern(l)));   %'MarkerFaceColor',rand(1,3));
% end

