function [F_BBu,W_BBu]=UPBB_general_BDSVD(F_RFu,W_RF,H_user,U,Ns_alloc,N_sub,Nrf_per_MS,Nrf_BS)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For digital domain beamfromer of hybrid beamformer
%   Apply block diagonalization(BD) to handle the inter-user interference and
%   SVD scheme to cancel inter-stream interference(within the same user)
%   Same data stream allocation for entire frequency tones
%--------------------------intput parameter--------------------------------
%   F_RFu       :Analog precoder at each user
%   W_RF        :Analog combiner at BS
%   H_user      :Uplink channel
%   U           :Number of user
%   Ns_alloc    :Number of data streams for each users
%   N_sub       :Number of subcarriers
%   Nrf_per_MS  :Number of RF chanis per user 
%   Nrf_BS      :Number of RF chanis at BS
%--------------------------output parameter--------------------------------
%   F_BBu       :Digital precoder at each user
%   W_BBu       :Digital combiner at BS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    W_BBu = zeros(Nrf_BS,max(Ns_alloc),U,N_sub);
    F_BBu = zeros(Nrf_per_MS,max(Ns_alloc),U,N_sub);
    Equbb_ch = zeros(Nrf_BS,Nrf_per_MS,U,N_sub);
    Fj = zeros(Nrf_per_MS,max(Ns_alloc),U,N_sub);
    rank_int_ch = zeros(U,N_sub);
    % Coordinated Tx-Rx BD Algorithm
    for sub_index=1:1:N_sub
        Hu = [];
        for u=1:1:U
            Equbb_ch(:,:,u,sub_index) = W_RF'*H_user(:,:,sub_index,u)*F_RFu(:,:,u);
            [~,~,guess_V] = svd(Equbb_ch(:,:,u,sub_index));
            Fj(:,1:Ns_alloc(u),u,sub_index) = guess_V(:,1:Ns_alloc(u));
            Hu = [Hu,Equbb_ch(:,:,u,sub_index)*Fj(:,1:Ns_alloc(u),u,sub_index)];
        end

        for u=1:1:U
            Hu_tilte = Hu;
            Hu_tilte(:,sum(Ns_alloc(1:1:u-1))+1:1:sum(Ns_alloc(1:1:u))) = [];

            [U_tilte,S_tilte,V_tilte] = svd(Hu_tilte);       
            rank_int_ch(u,sub_index) = rank(Hu_tilte);

            [UU_tilte,SS_tilte,VV_tilte] = svd(U_tilte(:,rank_int_ch(u,sub_index)+1:end)'*Equbb_ch(:,:,u,sub_index)*Fj(:,1:Ns_alloc(u),u,sub_index));
            F_BBu(:,1:Ns_alloc(u),u,sub_index) = Fj(:,1:Ns_alloc(u),u,sub_index)*VV_tilte(:,1:Ns_alloc(u));
            W_BBu(:,1:Ns_alloc(u),u,sub_index) = U_tilte(:,rank_int_ch(u,sub_index)+1:end)*UU_tilte(:,1:Ns_alloc(u));

        end
    end
    
end % end of function 