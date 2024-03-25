function [R,R_iku,S,inter_stream_I,inter_user_I,I,N]=Multi_Ns_UP_rate(P_uik,F_HBFu,W_HBFu,H_user,N0,U,Ns_alloc,N_sub)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   For hybrid beamformer
%   The ahievable rate for upink multi-ant and multi data stream transmission 
%--------------------------intput parameter--------------------------------
%   P_uik       :Transmit power of the i-th stream of the u-th user at the k-th subacarrier
%   F_HBFu      :Hybrid precoder
%   W_HBFu      :Hybrid cobiner
%   H_user      :Uplink channel
%   N0          :Noise power at k-th subacarrier
%   U           :Number of user
%   Ns_alloc    :Number of data streams for each users
%   N_sub       :Number of subcarriers
%--------------------------output parameter--------------------------------
%   R             :Ahievable rate for entire system(bits/sec/Hz)
%   R_iku         :Ahievable rate of the i-th stream of the u-th user at the k-th subacarrier(bits/sec/Hz)
%   S             :Desire signal power
%   inter_stream_I:Inter-stream interference(within the same user)
%   inter_user_I  :Inter-user interference
%   I             :Total interference(inter_stream_I+inter_user_I)
%   N             :Effective noise power
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    HBF_Rsum = 0;

    S = zeros(max(Ns_alloc),N_sub,U);
    inter_stream_I = zeros(max(Ns_alloc),N_sub,U);
    inter_user_I = zeros(max(Ns_alloc),N_sub,U);
    I = zeros(max(Ns_alloc),N_sub,U);
    N = zeros(max(Ns_alloc),N_sub,U);
    R_iku = zeros(max(Ns_alloc),N_sub,U);
    
    alluser_ind = 1:1:U;
    for sub_index=1:1:N_sub
          for u=1:1:U
              intuser_ind = alluser_ind;
              intuser_ind(:,u) = [];
              allds_ind = 1:1:Ns_alloc(u);
              
              for i=1:1:Ns_alloc(u)
                intds_ind = allds_ind;
                intds_ind(:,i) = [];

                S(i,sub_index,u) = P_uik(i,sub_index,u)*abs(W_HBFu(:,i,u,sub_index)'*H_user(:,:,sub_index,u)*F_HBFu(:,i,u,sub_index))^2;

                inter_stream_I(i,sub_index,u) = abs(W_HBFu(:,i,u,sub_index)'*H_user(:,:,sub_index,u)*F_HBFu(:,intds_ind,u,sub_index)).^2*P_uik(intds_ind,sub_index,u); 

                for m=intuser_ind
                    inter_user_I(i,sub_index,u) = inter_user_I(i,sub_index,u)+abs(W_HBFu(:,i,u,sub_index)'*H_user(:,:,sub_index,m)*F_HBFu(:,1:Ns_alloc(m),m,sub_index)).^2*P_uik(1:Ns_alloc(m),sub_index,m);
                end

                I(i,sub_index,u) = inter_stream_I(i,sub_index,u)+inter_user_I(i,sub_index,u);
                N(i,sub_index,u) = N0*norm(W_HBFu(:,i,u,sub_index))^2;

                R_iku(i,sub_index,u) = log2(1+S(i,sub_index,u)/(I(i,sub_index,u)+N(i,sub_index,u)));
                HBF_Rsum = HBF_Rsum+R_iku(i,sub_index,u);
              end % end of data stream index for u-th user                  
          end % end of user index 
    end % end of subcarrier index 
    R = HBF_Rsum/N_sub;
end