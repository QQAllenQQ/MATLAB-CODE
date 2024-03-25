function [ mu_kappa1, mu_kappa2 ] = functionComputeDistanceExpectation( pathLossExp_vec,distancesPathLossExp_vec,Gn_vec,lambda )

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
%This function returns the two auxiliary terms $\mu_1$ and $\mu_2$ used in
%Lemma 1 in the article to compute the UatF bound on the uplink average
%ergodic spectral efficiency.
%
%
%INPUT:
%pathLossExp_vec            = Path loss exponents used in the multi slope model
%distancesPathLossExp_vec   = Reference distances at which a change in power 
%                           decandence occur used in the multi slope model
%Gn_vec                     = multislope pathloss coefficients 
%lambdaRange                = Density of BSs for the random deployment.
%
%OUTPUT:
%mu_kappa1      = 1 x 1 auxiliary term in Lemma 1 (\mu_1)
%mu_kappa2      = 1 x 1 auxiliary term in Lemma 1 (\mu_2)

%
N = length(pathLossExp_vec);
alpha_vec = fliplr([pathLossExp_vec(1), pathLossExp_vec]);
Gn_vec = fliplr([Gn_vec(1), Gn_vec]);
R_vec = [distancesPathLossExp_vec./1e3 1e10];

%Compute the \mu term in Appendix B
mu_kappa1_sum = 0;
mu_kappa2_sum = 0;
%Go through all the number of slopes
for n = 1:N
    coeff1sum = 0;
    coeff2sum = 0;
    for ii = n+1:N+1
        Ri_1 = R_vec(ii);
        Ri = R_vec(ii+1);
        alphai = alpha_vec(ii);
        coeff1sum = coeff1sum + (Gn_vec(ii)/Gn_vec(n)) * ((Ri_1^(2-alphai)-Ri^(2-alphai))/(alphai-2));
        coeff2sum = coeff2sum + (Gn_vec(ii)/Gn_vec(n))^2 * ((Ri_1^(2-2*alphai)-Ri^(2-2*alphai))/(2*alphai-2));
    end
    Rn_1 = R_vec(n);
    Rn = R_vec(n+1);
    alphan = alpha_vec(n);
    cn_kappa1 = coeff1sum - (Rn^(2-alphan))/(alphan-2);
    cn_kappa2 = coeff2sum - (Rn^(2-2*alphan))/(2*alphan-2);
    %
    gamman_diff_kappa1 = ( gammainc(pi*lambda*Rn_1^2,1+alphan/2,'upper')-gammainc(pi*lambda*Rn^2,1+alphan/2,'upper') )*gamma(1+alphan/2);
    gamman_diff_kappa2 = ( gammainc(pi*lambda*Rn_1^2,1+alphan,'upper')-gammainc(pi*lambda*Rn^2,1+alphan,'upper') )*gamma(1+alphan);
    gamman_diff = ( gammainc(pi*lambda*Rn_1^2,2,'upper')-gammainc(pi*lambda*Rn^2,2,'upper') )*gamma(2);
    %
    mu_kappa1_sum = mu_kappa1_sum + ( 2*gamman_diff/(alphan-2) + (2*pi*lambda/(pi*lambda)^(alphan/2))*cn_kappa1*gamman_diff_kappa1 );
    mu_kappa2_sum = mu_kappa2_sum + ( 2*gamman_diff/(2*alphan-2) + (2*pi*lambda/(pi*lambda)^alphan)*cn_kappa2*gamman_diff_kappa2 );
end

%save
mu_kappa1 = mu_kappa1_sum;
mu_kappa2 = mu_kappa2_sum;



end

