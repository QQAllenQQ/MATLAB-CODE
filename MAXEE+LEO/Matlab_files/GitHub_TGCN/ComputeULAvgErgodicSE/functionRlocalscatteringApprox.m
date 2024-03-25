function R = functionRlocalscatteringApprox(M,theta,ASDdeg,antennaSpacing)

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
%This functions generates the spatial correlation matrix for the local scattering model
%with Gaussian angular distribution. The small-angle approximation described in Section 2.6.2 
%is used to increase efficiency, thus this function should only be used for ASDs below 15 degrees.
%For implementation details, see Section 2.6.2 in the monohraph:
%
% Emil Björnson, Jakob Hoydis and Luca Sanguinetti, "Massive MIMO Networks: Spectral, Energy, and Hardware Efficiency", 
% Foundations and Trends® in Signal Processing: Vol. 11: No. 3-4, pp 154-655. http://dx.doi.org/10.1561/2000000093
%
%
%INPUT:
%M              = Number of antennas
%theta          = Nominal angle
%ASDdeg         = Angular standard deviation around the nominal angle
%                 (measured in degrees)
%antennaSpacing = (Optional) Spacing between antennas (in wavelengths)
%
%OUTPUT:
%R              = M x M spatial correlation matrix


%Set the antenna spacing if not specified by input
if  nargin < 4
    
    %Half a wavelength distance
    antennaSpacing = 1/2;
    
end


%Compute the ASD in radians based on input
ASD = ASDdeg*pi/180; 


%The correlation matrix has a Toeplitz structure, so we only need to
%compute the first row of the matrix
firstRow = zeros(M,1);


%Go through all the columns of the first row
for column = 1:M
    
    %Distance from the first antenna
    distance = column-1;
    
    %Compute the approximated integral as in (2.24)
    firstRow(column) = exp(1i*2*pi*antennaSpacing*sin(theta)*distance)*exp(-ASD^2/2 * ( 2*pi*antennaSpacing*cos(theta)*distance )^2);
    
end

%Compute the spatial correlation matrix by utilizing the Toeplitz structure
R = toeplitz(firstRow);
