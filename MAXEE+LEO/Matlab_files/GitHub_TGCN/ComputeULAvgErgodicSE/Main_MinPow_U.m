% Upink multi-ant 
% P vs U
tic;clc;clear;close all;
%% parameter seting
% BS parameter
N_BS=128;
Nrf_BS=12;                                                                  % number of RF chanis at BS(This parameter determine the most number of receive data stream)
Total_U = [20 30 40 50 60 70];

% Mobile parameter
U=length(Total_U);                                                                   % number of users
Ns_per_MS=2;                                                               % number of data streams per user
Nrf_per_MS=3;                                                              % number of RF chanis per user(This parameter determine the most number of transmitt data stream per user) 
N_per_MS=16;

h_BS=10;                                                                   % BS height[meter]
R=100;                                                                     % distance of BS and user[meter]
d0=10;                                                                     % minimum  distance between BS and UE(meter)
fc=28;                                                                     % carrier frequency(28 GHz)
N_sub=16;                                                                  % number of subcarriers
N_bits = 75;

Ns_TU_alloc = Ns_per_MS*ones(1,U);
Ns_TU_alloc_3 = Ns_per_MS*ones(1,U);
Ns_FD_alloc = Ns_per_MS*ones(1,U);
Ns_CHAVG_alloc = Ns_per_MS*ones(1,U);
Ns_equal_alloc_2 = Ns_per_MS*ones(1,U);
Ns_TU_alloc_4 = Ns_per_MS*ones(1,U);
Ns_FD_alloc_2 = Ns_per_MS*ones(1,U);
Ns_CHAVG_alloc_2 = Ns_per_MS*ones(1,U);


System_BW=(20*10^6);                                                       % System Bandwidth(50 MHz)
sub_BW=System_BW/N_sub;                                                    % Subcarrier frequency spacing
noise_variance_dBm=-159+10*log10(sub_BW);                                  % noise density=-159 dBm/Hz
N0=10^(noise_variance_dBm/10)*10^-3;                                       % noise power[W] at each subcarrier

% Ref:5G;NR; User Equipment (UE) radio transmission and reception;Part 2:
% Range 2 Standalone (3GPP TS 38.101-2 version 15.3.0 Release 15) CH 6.
UE_Max_Tx_dBm=30;                                                          % UE maximum output power limits for power class 3(dBm)
P_Max=(10^(UE_Max_Tx_dBm/10)*10^-3);                                       % UE maximum output power limits for power class 3(W) on entire bandwidth
UE_TX_POW_dBm=20;
ch_reali=10^8;   
P_total_user=(10.^(UE_TX_POW_dBm/10)/1000);                                                                                                                       
ang_spread=10*pi/180;
Nray=10;
PT_final1 = zeros(ch_reali,U);
PT_final1_3 = zeros(ch_reali,U);
PT_final1_4 = zeros(ch_reali,U);
PT_final2 = zeros(ch_reali,U);
PT_final3 = zeros(ch_reali,U);                                                                      
PT_final4 = zeros(ch_reali,U);
PT_final5 = zeros(ch_reali,U);                                                                      
PT_final6 = zeros(ch_reali,U);
%% MODCOD for DVB-S2 (Switching level v.s. SINR)
% Modulation x code rate
data_rate_st = [0    0.4902   1.1883   2.4786   3.5673   4.4530];                       %     X   MCS1   MCS2   MCS3   MCS4   MCS5
data_rate_ed = [0.4902  1.1883   2.4786  3.5673   4.4530    Inf];                       %  MCS1   MCS2   MCS3   MCS4   MCS5

% SINR threshold for each MCS index at 10% BLER
switching_level_st_dB = [-Inf  -2.54  2.07  9.28  12.98  15.76];             %     X   MCS1   MCS2   MCS3   MCS4   MCS5
switching_level_ed_dB = [-2.54  2.07  9.28  12.98  15.76 Inf];             %  MCS1   MCS2   MCS3   MCS4   MCS5
switching_level_st = 10.^(switching_level_st_dB/10);
switching_level_ed = 10.^(switching_level_ed_dB/10);
 
data_rate_change = data_rate_ed-data_rate_st;                                              % Date rate increased    [0.3770    0.7988    1.2305    1.4960    1.6524    Inf]
switching_level_change = switching_level_ed-switching_level_st;                            % Switch increased power [0.3665    1.2175    4.8120   15.9361   61.1416    Inf]
%% main program
for ch=1:1:ch_reali
    for r=1:1:U 
        [ch r]
        if Nrf_BS<=Nrf_per_MS*U
            Ns_total = Nrf_BS;
        else
            Ns_total = Nrf_per_MS*U;
        end   
              %% UE are randomly dropped in the cell over distances ranging between 10 m and the cell radius
        user_location_phase=pi*(rand(1,Total_U(r))*2-1);%-pi~pi
        user_location_radius=rand(1,Total_U(r))*(R-d0)+d0; 
        user_location=user_location_radius.*exp(j*user_location_phase);
        x_location=real(user_location);
        y_location=imag(user_location);
              %% clustered channel model
        H_user = zeros(N_BS,N_per_MS,N_sub,Total_U(r));
        for uu=1:1:Total_U(r)
            d_2D(uu)=norm([x_location(uu) y_location(uu)],2);                     % unit is meter    
            [PL(uu),path_loss_dB(uu)]=PathLoss(fc,h_BS,d_2D(uu));                            % path loss model
            [CH] =Satellite_channel(N_per_MS,N_BS,N_sub,1,Nray,ang_spread,path_loss_dB(uu),fc);
            H_user(:,:,:,uu)=CH(:,:,:,1)/sqrt(PL(uu));
        end
        %% Scheme 1:  Proposed User Selection
           max_userNs_1=4;
           P = N_BS;
           max_singular_1 = zeros(N_sub,Total_U(U),1);
        for nuu=1:Total_U(r)
           for ncc = 1:N_sub
              [~,sl_1(:,:,ncc),~] = svd(H_user(:,:,ncc,nuu),'econ'); 
              sl_(:,:,ncc) = diag(sl_1(:,:,ncc));               
           end
           for chh = 1:1
               max_singular_1(:,nuu,chh) = reshape(sl_(chh,1,:),1,N_sub);
               [sorted_gain(:,:,chh),IX(:,:,chh)] = sort(max_singular_1(:,:,chh),2,'descend');
           end          
        end

        selected_subcarriers = IX(:,1:U);

% find maximum channel        
if max_userNs_1 == 1
    for nc = 1:N_sub
        for nuu = 1:Total_U(r)
            if selected_subcarriers(nc,nuu)~=0
                channel_max(nc,1) = selected_subcarriers(nc,nuu);
                break;
            end
        end
    end
else
    %---------- submax channel ------------%
    for nsc = 1:N_sub
        i=1;
        for nu = 1:Total_U(r)
            if selected_subcarriers(nsc,nu)~=0
                channel_max(nsc,i) = selected_subcarriers(nsc,nu);
                i=i+1;
            end
            if i == max_userNs_1+1       % Ns=2 - i=2
                break;
            end
        end
    end 
end
selected_subcarriers = sort(selected_subcarriers,2,'descend');    
Sel_User_ind = selected_subcarriers(1,1:U);
       %% 這邊是為了把被選到的使用者通道挑出來存放
    Sel_H_user = zeros(N_BS,N_per_MS,N_sub,U);
    for u=1:1:U
        Sel_H_user(:,:,:,u) = H_user(:,:,:,Sel_User_ind(u));     
    end
    %% TUMD
        % analog precoder at MS via Tensor unfolding
        long_unfold_Hu = zeros((N_sub)*N_BS,N_per_MS,U); 
        long_unfold_F_RFu = zeros(N_per_MS,Nrf_per_MS,U);
        for u=1:1:U
            for sub_index=1:1:N_sub
                long_unfold_Hu((sub_index-1)*N_BS+1:1:sub_index*N_BS,:,u) = Sel_H_user(:,:,sub_index,u);
            end
            [long_unfold_V,long_unfold_D] = eig(long_unfold_Hu(:,:,u)'*long_unfold_Hu(:,:,u));  
            [~,longunf_ind] = sort(diag(long_unfold_D),'descend');
            long_unfold_F_RFu(:,:,u) = extract_phase_v2(long_unfold_V(:,longunf_ind(1:1:Nrf_per_MS)),N_per_MS);
        end    

        % analog combiner at BS via Tensor unfolding and RF chain assignment
        horiz_unfold_Hu = zeros(N_BS,Nrf_per_MS*(N_sub),U);
        horiz_unfold_W_RF = [];
        eigval_buff = zeros(N_BS,U);
        eigvec_buff = zeros(N_BS,N_BS,U);

        for u=1:1:U
            for sub_index=1:1:N_sub
                horiz_unfold_Hu(:,(sub_index-1)*Nrf_per_MS+1:1:sub_index*Nrf_per_MS,u) = Sel_H_user(:,:,sub_index,u)*long_unfold_F_RFu(:,:,u);
            end
            [horiz_unfold_V,horiz_unfold_D] = eig(horiz_unfold_Hu(:,:,u)*horiz_unfold_Hu(:,:,u)');
            [horizunfout,horizunfind] = sort(diag(horiz_unfold_D),'descend');

            eigval_buff(:,u) = real(horizunfout);
            eigvec_buff(:,:,u) = horiz_unfold_V(:,horizunfind);  
        end           
        for u=1:1:U
            horiz_unfold_W_RF = [horiz_unfold_W_RF,extract_phase_v2(eigvec_buff(:,1:1:Ns_TU_alloc(u),u),N_BS)];
        end   
        % BD(ignore the noise)+SVD beamforming based on effective channel
        [TU_F_BBu,TU_W_BBu]=UPBB_general_BDSVD(long_unfold_F_RFu,horiz_unfold_W_RF,Sel_H_user,U,Ns_TU_alloc,N_sub,Nrf_per_MS,Nrf_BS);

        % Compute hybrid beamforming matrix
        TU_F_HBFu = zeros(N_per_MS,max(Ns_TU_alloc),U,N_sub);
        TU_W_HBFu = zeros(N_BS,max(Ns_TU_alloc),U,N_sub);
        for sub_index=1:1:N_sub       
            for u=1:1:U
                TU_F_HBFu(:,1:Ns_TU_alloc(u),u,sub_index) = long_unfold_F_RFu(:,:,u)*TU_F_BBu(:,1:Ns_TU_alloc(u),u,sub_index);

                % normalize transmit HBF Weight vector to one
                for ns = 1:1:Ns_TU_alloc(u)
                    TU_F_HBFu(:,ns,u,sub_index) = TU_F_HBFu(:,ns,u,sub_index)/norm(TU_F_HBFu(:,ns,u,sub_index)); 
                end

                TU_W_HBFu(:,1:Ns_TU_alloc(u),u,sub_index) = horiz_unfold_W_RF*TU_W_BBu(:,1:Ns_TU_alloc(u),u,sub_index);
            end
        end
        %% Effective channel calculation
    TU_H_eff = zeros(N_BS,Ns_total,N_sub);
    TU_W_HBF_fixed = zeros(N_BS,Ns_total,N_sub);
    TU_result_CH = zeros(Ns_total,N_sub);
    Noise_eff = zeros(Ns_total,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            for nn = 1:Ns_TU_alloc(u)  
            TU_W_HBF_fixed(:,sum(Ns_TU_alloc(1:1:u-1))+1:1:sum(Ns_TU_alloc(1:1:u)),sub_index) = TU_W_HBFu(:,1:Ns_TU_alloc(u),u,sub_index);
            TU_H_eff(:,sum(Ns_TU_alloc(1:1:u-1))+1:1:sum(Ns_TU_alloc(1:1:u)),sub_index) = Sel_H_user(:,:,sub_index,u)*TU_F_HBFu(:,1:Ns_TU_alloc(u),u,sub_index);
            
            Noise_eff(sum(Ns_TU_alloc(1:1:u-1))+nn,sub_index) = N0*norm(TU_W_HBF_fixed(sum(Ns_TU_alloc(1:1:u-1))+nn,sum(Ns_TU_alloc(1:1:u-1))+nn,sub_index)')^2;
            Noise_eff;
            TU_result_CH(sum(Ns_TU_alloc(1:1:u-1))+nn,sub_index) = norm(TU_W_HBF_fixed(sum(Ns_TU_alloc(1:1:u-1))+nn,sum(Ns_TU_alloc(1:1:u-1))+nn,sub_index)'*TU_H_eff(sum(Ns_TU_alloc(1:1:u-1))+nn,sum(Ns_TU_alloc(1:1:u-1))+nn,sub_index))^2;  %% Heff
            end
        end  
    end
    TU_result_CH;
         %% Step3. MODCOD index loading 
        Unit_bits = data_rate_change(1)*ones(Ns_total,N_sub);
        SINR = switching_level_change(1)*ones(Ns_total,N_sub); 
        MODCOD_Index_table1 = zeros(Ns_total,N_sub);                                                         % Set up MODCOD index table
        PT_unit_table1 = (SINR./Unit_bits.*Noise_eff)./TU_result_CH;   
        for nu = 1:U 
            st = sum(Ns_TU_alloc(1:(nu-1)))+1;
            ed = sum(Ns_TU_alloc(1:nu)); 
            bits = N_bits;                                                                             % Available bits                           
            while bits > 0
                [x,y] = find(PT_unit_table1(st:ed,:) == min(min(PT_unit_table1(st:ed,:))));                          % Find least power increase 
                MODCOD_Index_table1(st+x(1)-1,y(1)) = MODCOD_Index_table1(st+x(1)-1,y(1)) + 1;                             % Update MODCOD index table
                PT_unit_table1(st+x(1)-1,y(1)) = (switching_level_change(MODCOD_Index_table1(st+x(1)-1,y(1))+1)/data_rate_change(MODCOD_Index_table1(st+x(1)-1,y(1))+1)*Noise_eff(st+x(1)-1,y(1)))./TU_result_CH(st+x(1)-1,y(1));
                bits = bits - data_rate_change(MODCOD_Index_table1(st+x(1)-1,y(1)));                                    % Update transmitted power efficiency table                                                                                                                                                   % Update available bits
            end
        end 
        MODCOD_Index_table1;
        bits_table1 = data_rate_st(MODCOD_Index_table1+1);
        %% Calculate transmitted power
        for k = 1:N_sub
            for nf = 1:Ns_total
                PT_final_table1(nf,k) = (switching_level_st(MODCOD_Index_table1(nf,k)+1)*Noise_eff(nf,k))/TU_result_CH(nf,k);  
            end
        end 
        PT_final_table1;
        %% TUMD (K/4)
    % analog precoder at MS via Tensor unfolding
    long_unfold_Hu_3 = zeros((N_sub/4)*N_BS,N_per_MS,U); 
    long_unfold_F_RFu_3 = zeros(N_per_MS,Nrf_per_MS,U);
    for u=1:1:U
        for sub_index2=2:4:14
            long_unfold_Hu_3((sub_index2-1)*N_BS+1:1:sub_index2*N_BS,:,u) = Sel_H_user(:,:,sub_index2,u);
        end
        [long_unfold_V_3,long_unfold_D_3] = eig(long_unfold_Hu_3(:,:,u)'*long_unfold_Hu_3(:,:,u));  
        [~,longunf_ind_3] = sort(diag(long_unfold_D_3),'descend');
        long_unfold_F_RFu_3(:,:,u) = extract_phase_v2(long_unfold_V_3(:,longunf_ind_3(1:1:Nrf_per_MS)),N_per_MS);
    end    
    
    % analog combiner at BS via Tensor unfolding and RF chain assignment
    horiz_unfold_Hu_3 = zeros(N_BS,Nrf_per_MS*(N_sub/4),U);
    horiz_unfold_W_RF_3 = [];
    eigval_buff_3 = zeros(N_BS,U);
    eigvec_buff_3 = zeros(N_BS,N_BS,U);
    
    for u=1:1:U
        for sub_index2=2:4:14
            horiz_unfold_Hu_3(:,(sub_index2-1)*Nrf_per_MS+1:1:sub_index2*Nrf_per_MS,u) = Sel_H_user(:,:,sub_index2,u)*long_unfold_F_RFu_3(:,:,u);
        end
        [horiz_unfold_V_3,horiz_unfold_D_3] = eig(horiz_unfold_Hu_3(:,:,u)*horiz_unfold_Hu_3(:,:,u)');
        [horizunfout_3,horizunfind_3] = sort(diag(horiz_unfold_D_3),'descend');
        
        eigval_buff_3(:,u) = real(horizunfout_3);
        eigvec_buff_3(:,:,u) = horiz_unfold_V_3(:,horizunfind_3);  
    end   
 
    for u=1:1:U
        horiz_unfold_W_RF_3 = [horiz_unfold_W_RF_3,extract_phase_v2(eigvec_buff_3(:,1:1:Ns_TU_alloc_3(u),u),N_BS)];
    end
    
    % BD(ignore the noise)+SVD beamforming based on effective channel
    [TU_F_BBu_3,TU_W_BBu_3]=UPBB_general_BDSVD(long_unfold_F_RFu_3,horiz_unfold_W_RF_3,Sel_H_user,U,Ns_TU_alloc_3,N_sub,Nrf_per_MS,Nrf_BS);
    
    % Compute hybrid beamforming matrix
    TU_F_HBFu_3 = zeros(N_per_MS,max(Ns_TU_alloc_3),U,N_sub);
    TU_W_HBFu_3 = zeros(N_BS,max(Ns_TU_alloc_3),U,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            TU_F_HBFu_3(:,1:Ns_TU_alloc_3(u),u,sub_index) = long_unfold_F_RFu_3(:,:,u)*TU_F_BBu_3(:,1:Ns_TU_alloc_3(u),u,sub_index);
            
            % normalize transmit HBF Weight vector to one
            for ns = 1:1:Ns_TU_alloc_3(u)
                TU_F_HBFu_3(:,ns,u,sub_index) = TU_F_HBFu_3(:,ns,u,sub_index)/norm(TU_F_HBFu_3(:,ns,u,sub_index)); 
            end
            
            TU_W_HBFu_3(:,1:Ns_TU_alloc_3(u),u,sub_index) = horiz_unfold_W_RF_3*TU_W_BBu_3(:,1:Ns_TU_alloc_3(u),u,sub_index);
        end
    end
    %% Effective channel calculation
    TU_H_eff_3 = zeros(N_BS,Ns_total,N_sub);
    TU_W_HBF_fixed_3 = zeros(N_BS,Ns_total,N_sub);
    TU_result_CH_3 = zeros(Ns_total,N_sub);
    Noise_eff_3 = zeros(Ns_total,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            for nn_3 = 1:Ns_TU_alloc_3(u)
            TU_W_HBF_fixed_3(:,sum(Ns_TU_alloc_3(1:1:u-1))+1:1:sum(Ns_TU_alloc_3(1:1:u)),sub_index) = TU_W_HBFu_3(:,1:Ns_TU_alloc_3(u),u,sub_index);
            TU_H_eff_3(:,sum(Ns_TU_alloc_3(1:1:u-1))+1:1:sum(Ns_TU_alloc_3(1:1:u)),sub_index) = Sel_H_user(:,:,sub_index,u)*TU_F_HBFu_3(:,1:Ns_TU_alloc_3(u),u,sub_index);
            Noise_eff_3(sum(Ns_TU_alloc_3(1:1:u-1))+nn_3,sub_index) = N0*norm(TU_W_HBF_fixed_3(sum(Ns_TU_alloc_3(1:1:u-1))+nn_3,sum(Ns_TU_alloc_3(1:1:u-1))+nn_3,sub_index)')^2;
            Noise_eff_3;
            TU_result_CH_3(sum(Ns_TU_alloc_3(1:1:u-1))+nn_3,sub_index) = norm(TU_W_HBF_fixed_3(sum(Ns_TU_alloc_3(1:1:u-1))+nn_3,sum(Ns_TU_alloc_3(1:1:u-1))+nn_3,sub_index)'*TU_H_eff_3(sum(Ns_TU_alloc_3(1:1:u-1))+nn_3,sum(Ns_TU_alloc_3(1:1:u-1))+nn_3,sub_index))^2;  %% Heff
            end  
        end
    end
    TU_result_CH_3;
         %% Step3. MODCOD index loading 
        Unit_bits = data_rate_change(1)*ones(Ns_total,N_sub);
        SINR = switching_level_change(1)*ones(Ns_total,N_sub); 
        MODCOD_Index_table1_3 = zeros(Ns_total,N_sub);                                                         % Set up MODCOD index table
        PT_unit_table1_3 = (SINR./Unit_bits.*Noise_eff_3)./TU_result_CH_3;   
        for nu = 1:U 
            st = sum(Ns_TU_alloc_3(1:(nu-1)))+1;
            ed = sum(Ns_TU_alloc_3(1:nu)); 
            bits = N_bits;                                                                             % Available bits                           
            while bits > 0
                [x,y] = find(PT_unit_table1_3(st:ed,:) == min(min(PT_unit_table1_3(st:ed,:))));                          % Find least power increase 
                MODCOD_Index_table1_3(st+x(1)-1,y(1)) = MODCOD_Index_table1_3(st+x(1)-1,y(1)) + 1;                             % Update MODCOD index table
                PT_unit_table1_3(st+x(1)-1,y(1)) = (switching_level_change(MODCOD_Index_table1_3(st+x(1)-1,y(1))+1)/data_rate_change(MODCOD_Index_table1_3(st+x(1)-1,y(1))+1)*Noise_eff_3(st+x(1)-1,y(1)))./TU_result_CH_3(st+x(1)-1,y(1));
                bits = bits - data_rate_change(MODCOD_Index_table1_3(st+x(1)-1,y(1)));                                    % Update transmitted power efficiency table                                                                                                                                                   % Update available bits
            end
        end 
        MODCOD_Index_table1_3;
        bits_table1_3 = data_rate_st(MODCOD_Index_table1_3+1);
        %% Calculate transmitted power
        for k = 1:N_sub
            for nf = 1:Ns_total
                PT_final_table1_3(nf,k) = (switching_level_st(MODCOD_Index_table1_3(nf,k)+1)*Noise_eff_3(nf,k))/TU_result_CH_3(nf,k);  
            end
        end 
        PT_final_table1_3;
       %% Fully-digital beamforming
        FD_Hu = zeros(N_BS,N_per_MS*(U),N_sub); 
        FD_U_tilte = zeros(N_BS,N_BS,U,N_sub);
        rank_int_ch = zeros(U,N_sub);
        for sub_index=1:1:N_sub
            for u=1:1:U
                FD_Hu(:,(u-1)*N_per_MS+1:1:u*N_per_MS,sub_index) = Sel_H_user(:,:,sub_index,u);
            end
            for u=1:1:U
                FD_Hu_tilte = FD_Hu(:,:,sub_index);
                FD_Hu_tilte(:,(u-1)*N_per_MS+1:1:u*N_per_MS) = [];

                [FD_U_tilte(:,:,u,sub_index),~,~] = svd(FD_Hu_tilte);      
                rank_int_ch(u,sub_index) = rank(FD_Hu_tilte);
            end     
        end
        FD_W_u = zeros(N_BS,max(Ns_FD_alloc),U,N_sub);
        FD_F_u = zeros(N_per_MS,max(Ns_FD_alloc),U,N_sub);
        for sub_index=1:1:N_sub
            for u=1:1:U
                [FD_UU_tilte,~,FD_VV_tilte] = svd(FD_U_tilte(:,rank_int_ch(u,sub_index)+1:end,u,sub_index)'*Sel_H_user(:,:,sub_index,u));
                FD_F_u(:,1:Ns_FD_alloc(u),u,sub_index) = FD_VV_tilte(:,1:Ns_FD_alloc(u));
                FD_W_u(:,1:Ns_FD_alloc(u),u,sub_index) = FD_U_tilte(:,rank_int_ch(u,sub_index)+1:end,u,sub_index)*FD_UU_tilte(:,1:Ns_FD_alloc(u));    
            end     
        end
        %% Effective channel calculation
    FD_H_eff = zeros(N_BS,Ns_total,N_sub);
    FD_W_HBF_fixed = zeros(N_BS,Ns_total,N_sub);
    FD_result_CH = zeros(Ns_total,N_sub);
    FD_Noise_eff = zeros(Ns_total,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            for fn = 1:Ns_FD_alloc(u)
            FD_W_HBF_fixed(:,sum(Ns_FD_alloc(1:1:u-1))+1:1:sum(Ns_FD_alloc(1:1:u)),sub_index) = FD_W_u(:,1:Ns_FD_alloc(u),u,sub_index);
            FD_H_eff(:,sum(Ns_FD_alloc(1:1:u-1))+1:1:sum(Ns_FD_alloc(1:1:u)),sub_index) = Sel_H_user(:,:,sub_index,u)*FD_F_u(:,1:Ns_FD_alloc(u),u,sub_index);
            
            FD_Noise_eff(sum(Ns_FD_alloc(1:1:u-1))+fn,sub_index) = N0*norm(FD_W_HBF_fixed(sum(Ns_FD_alloc(1:1:u-1))+fn,sum(Ns_FD_alloc(1:1:u-1))+fn,sub_index)')^2;
            FD_result_CH(sum(Ns_FD_alloc(1:1:u-1))+fn,sub_index) = norm(FD_W_HBF_fixed(sum(Ns_FD_alloc(1:1:u-1))+fn,sum(Ns_FD_alloc(1:1:u-1))+fn,sub_index)'*FD_H_eff(sum(Ns_FD_alloc(1:1:u-1))+fn,sum(Ns_FD_alloc(1:1:u-1))+fn,sub_index))^2;  %% Heff
            end
        end
    end
    FD_result_CH;
         %% Step3. MODCOD index loading 
        Unit_bits = data_rate_change(1)*ones(Ns_total,N_sub);
        SINR = switching_level_change(1)*ones(Ns_total,N_sub); 
        MODCOD_Index_table2 = zeros(Ns_total,N_sub);                                                         % Set up MODCOD index table
        PT_unit_table2 = (SINR./Unit_bits.*FD_Noise_eff)./FD_result_CH;                 % Set up transmitted power efficiency table    
        for nu = 1:U
            st = sum(Ns_FD_alloc(1:1:nu-1))+1;
            ed = sum(Ns_FD_alloc(1:1:nu)); 
            bits = N_bits;                                                                             % Available bits                           
            while bits > 0
                [x,y] = find(PT_unit_table2(st:ed,:) == min(min(PT_unit_table2(st:ed,:))));                          % Find least power increase 
                MODCOD_Index_table2(st+x(1)-1,y(1)) = MODCOD_Index_table2(st+x(1)-1,y(1)) + 1;                             % Update MCS index table
                PT_unit_table2(st+x(1)-1,y(1)) = (switching_level_change(MODCOD_Index_table2(st+x(1)-1,y(1))+1)/data_rate_change(MODCOD_Index_table2(st+x(1)-1,y(1))+1)*FD_Noise_eff(st+x(1)-1,y(1)))./FD_result_CH(st+x(1)-1,y(1));
                bits = bits - data_rate_change(MODCOD_Index_table2(st+x(1)-1,y(1)));                                    % Update transmitted power efficiency table                                                                                                                                                   % Update available bits
            end
        end 
        MODCOD_Index_table2;
        bits_table2 = data_rate_st(MODCOD_Index_table2+1);
        %% Calculate transmitted power
        for k = 1:N_sub
            for nf1 = 1:Ns_total
                PT_final_table2(nf1,k) = (switching_level_st(MODCOD_Index_table2(nf1,k)+1)*FD_Noise_eff(nf1,k))/FD_result_CH(nf1,k);
            end
        end 
        PT_final_table2;
        %% Benchmark Scheme--ACMD
        % analog precoder at MS and analog combiner at BS via channel average
        chavg_F_RFu = zeros(N_per_MS,Nrf_per_MS,U);
        chavg_W_RF = [];

        chavg_sinval_buff = zeros(N_per_MS,U);
        chavg_L_sinvec_buff = zeros(N_BS,N_per_MS,U);
        % analog precoder
        for u=1:1:U
            Q = 0;
            for sub_index=1:1:N_sub
                Q = Q+Sel_H_user(:,:,sub_index,u)/N_sub;
            end
            [chavg_U,~,chavg_V] = svd(Q);
            chavg_F_RFu(:,:,u) = extract_phase_v2(chavg_V(:,1:Nrf_per_MS),N_per_MS);
            chavg_sinval_buff(:,u) = svd(Q);
            chavg_L_sinvec_buff(:,:,u) = chavg_U(:,1:N_per_MS);
        end  
        for u=1:1:U
            chavg_W_RF = [chavg_W_RF,extract_phase_v2(chavg_L_sinvec_buff(:,1:1:Ns_CHAVG_alloc(u),u),N_BS)];
        end

        % BD(ignore the noise)+SVD beamforming based on effective channel
        [chavg_F_BBu,chavg_W_BBu]=UPBB_general_BDSVD(chavg_F_RFu,chavg_W_RF,Sel_H_user,U,Ns_CHAVG_alloc,N_sub,Nrf_per_MS,Nrf_BS);

        % Compute hybrid beamforming matrix
        chavg_F_HBFu = zeros(N_per_MS,max(Ns_CHAVG_alloc),U,N_sub);
        chavg_W_HBFu = zeros(N_BS,max(Ns_CHAVG_alloc),U,N_sub);
        for sub_index=1:1:N_sub       
            for u=1:1:U
                chavg_F_HBFu(:,1:Ns_CHAVG_alloc(u),u,sub_index) = chavg_F_RFu(:,:,u)*chavg_F_BBu(:,1:Ns_CHAVG_alloc(u),u,sub_index);

                % normalize transmit HBF Weight vector to one
                for ns = 1:1:Ns_CHAVG_alloc(u)
                    chavg_F_HBFu(:,ns,u,sub_index) = chavg_F_HBFu(:,ns,u,sub_index)/norm(chavg_F_HBFu(:,ns,u,sub_index)); 
                end

                chavg_W_HBFu(:,1:Ns_CHAVG_alloc(u),u,sub_index) = chavg_W_RF*chavg_W_BBu(:,1:Ns_CHAVG_alloc(u),u,sub_index);
            end
        end
        %% Effective channel calculation
    chavg_H_eff = zeros(N_BS,Ns_total,N_sub);
    chavg_W_HBF_fixed = zeros(N_BS,Ns_total,N_sub);
    chavg_result_CH = zeros(Ns_total,N_sub);
    chavg_Noise_eff = zeros(Ns_total,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            for cn = 1:Ns_CHAVG_alloc(u)
            chavg_W_HBF_fixed(:,sum(Ns_CHAVG_alloc(1:1:u-1))+1:1:sum(Ns_CHAVG_alloc(1:1:u)),sub_index) = chavg_W_HBFu(:,1:Ns_CHAVG_alloc(u),u,sub_index);
            chavg_H_eff(:,sum(Ns_CHAVG_alloc(1:1:u-1))+1:1:sum(Ns_CHAVG_alloc(1:1:u)),sub_index) = Sel_H_user(:,:,sub_index,u)*chavg_F_HBFu(:,1:Ns_CHAVG_alloc(u),u,sub_index);
            
            chavg_Noise_eff(sum(Ns_CHAVG_alloc(1:1:u-1))+cn,sub_index) = N0*norm(chavg_W_HBF_fixed(sum(Ns_CHAVG_alloc(1:1:u-1))+cn,sum(Ns_CHAVG_alloc(1:1:u-1))+cn,sub_index)')^2;
            chavg_result_CH(sum(Ns_CHAVG_alloc(1:1:u-1))+cn,sub_index) = norm(chavg_W_HBF_fixed(sum(Ns_CHAVG_alloc(1:1:u-1))+cn,sum(Ns_CHAVG_alloc(1:1:u-1))+cn,sub_index)'*chavg_H_eff(sum(Ns_CHAVG_alloc(1:1:u-1))+cn,sum(Ns_CHAVG_alloc(1:1:u-1))+cn,sub_index))^2;  %% Heff
            end
        end
    end
    chavg_result_CH;
         %% Step3. MODCOD index loading 
        Unit_bits = data_rate_change(1)*ones(Ns_total,N_sub);
        SINR = switching_level_change(1)*ones(Ns_total,N_sub); 
        MODCOD_Index_table5 = zeros(Ns_total,N_sub);                                                         % Set up MODCOD index table
        PT_unit_table5 = (SINR./Unit_bits.*chavg_Noise_eff)./chavg_result_CH;                 % Set up transmitted power efficiency table     
        for nu = 1:U
            st = sum(Ns_CHAVG_alloc(1:1:nu-1))+1;
            ed = sum(Ns_CHAVG_alloc(1:1:nu)); 
            bits = N_bits;                                                                             % Available bits                           
            while bits > 0
                [x,y] = find(PT_unit_table5(st:ed,:) == min(min(PT_unit_table5(st:ed,:))));                          % Find least power increase 
                MODCOD_Index_table5(st+x(1)-1,y(1)) = MODCOD_Index_table5(st+x(1)-1,y(1)) + 1;                             % Update MODCOD index table
                PT_unit_table5(st+x(1)-1,y(1)) = (switching_level_change(MODCOD_Index_table5(st+x(1)-1,y(1))+1)/data_rate_change(MODCOD_Index_table5(st+x(1)-1,y(1))+1)*chavg_Noise_eff(st+x(1)-1,y(1)))./chavg_result_CH(st+x(1)-1,y(1));
                bits = bits - data_rate_change(MODCOD_Index_table5(st+x(1)-1,y(1)));                                    % Update transmitted power efficiency table                                                                                                                                                   % Update available bits
            end
        end 
        MODCOD_Index_table5;
        bits_table5 = data_rate_st(MODCOD_Index_table5+1);
        %% Calculate transmitted power
        for k = 1:N_sub
            for nf1 = 1:Ns_total
                PT_final_table5(nf1,k) = (switching_level_st(MODCOD_Index_table5(nf1,k)+1)*chavg_Noise_eff(nf1,k))/chavg_result_CH(nf1,k);
            end
        end 
        PT_final_table5;
        %% Scheme 2: Random User Selection
        % 從 使用者1 ~ 使用者Total_U 中選U個使用者來服務，Sel_User_ind輸出被服務使用者的index
        Sel_User_ind_2 = randperm(Total_U(r),U);
       %% 這邊是為了把被選到的使用者通道挑出來存放
    Sel_H_user_2 = zeros(N_BS,N_per_MS,N_sub,U);
    for u=1:1:U
        Sel_H_user_2(:,:,:,u) = H_user(:,:,:,Sel_User_ind_2(u));     
    end
    %% TUMD
        % analog precoder at MS via Tensor unfolding
        long_unfold_Hu_2 = zeros((N_sub)*N_BS,N_per_MS,U); 
        long_unfold_F_RFu_2 = zeros(N_per_MS,Nrf_per_MS,U);
        for u=1:1:U
            for sub_index=1:1:N_sub
                long_unfold_Hu_2((sub_index-1)*N_BS+1:1:sub_index*N_BS,:,u) = Sel_H_user_2(:,:,sub_index,u);
            end
            [long_unfold_V_2,long_unfold_D_2] = eig(long_unfold_Hu_2(:,:,u)'*long_unfold_Hu_2(:,:,u));  
            [~,longunf_ind_2] = sort(diag(long_unfold_D_2),'descend');
            long_unfold_F_RFu_2(:,:,u) = extract_phase_v2(long_unfold_V_2(:,longunf_ind_2(1:1:Nrf_per_MS)),N_per_MS);
        end    

        % analog combiner at BS via Tensor unfolding
        horiz_unfold_Hu_2 = zeros(N_BS,Nrf_per_MS*(N_sub),U);
        eigval_buff_2 = zeros(N_BS,U);
        eigvec_buff_2 = zeros(N_BS,N_BS,U);
    for u=1:1:U
        for sub_index=1:1:N_sub
            horiz_unfold_Hu_2(:,(sub_index-1)*Nrf_per_MS+1:1:sub_index*Nrf_per_MS,u) = Sel_H_user_2(:,:,sub_index,u)*long_unfold_F_RFu_2(:,:,u);
        end
        [horiz_unfold_V_2,horiz_unfold_D_2] = eig(horiz_unfold_Hu_2(:,:,u)*horiz_unfold_Hu_2(:,:,u)');
        [horizunfout_2,horizunfind_2] = sort(diag(horiz_unfold_D_2),'descend');
        
        eigval_buff_2(:,u) = real(horizunfout_2);
        eigvec_buff_2(:,:,u) = horiz_unfold_V_2(:,horizunfind_2);  
    end
    h_unf_W_RF_fixed_2 = [];
    for u=1:1:U
        h_unf_W_RF_fixed_2 = [h_unf_W_RF_fixed_2,extract_phase_v2(eigvec_buff_2(:,1:1:Ns_equal_alloc_2(u),u),N_BS)];
    end
    % BD(ignore the noise)+SVD beamforming based on effective channel
    [TU_F_BBu_fixed_2,TU_W_BBu_fixed_2]=UPBB_general_BDSVD(long_unfold_F_RFu_2,h_unf_W_RF_fixed_2,Sel_H_user_2,U,Ns_equal_alloc_2,N_sub,Nrf_per_MS,Nrf_BS);
    
    % Compute hybrid beamforming matrix
    TU_F_HBFu_fixed_2 = zeros(N_per_MS,max(Ns_equal_alloc_2),U,N_sub);
    TU_W_HBFu_fixed_2 = zeros(N_BS,max(Ns_equal_alloc_2),U,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            TU_F_HBFu_fixed_2(:,1:Ns_equal_alloc_2(u),u,sub_index) = long_unfold_F_RFu_2(:,:,u)*TU_F_BBu_fixed_2(:,1:Ns_equal_alloc_2(u),u,sub_index);
            % normalize transmit HBF Weight vector to one
            for ns1 = 1:1:Ns_equal_alloc_2(u)
                TU_F_HBFu_fixed_2(:,ns1,u,sub_index) = TU_F_HBFu_fixed_2(:,ns1,u,sub_index)/norm(TU_F_HBFu_fixed_2(:,ns1,u,sub_index)); 
            end
            
            TU_W_HBFu_fixed_2(:,1:Ns_equal_alloc_2(u),u,sub_index) = h_unf_W_RF_fixed_2*TU_W_BBu_fixed_2(:,1:Ns_equal_alloc_2(u),u,sub_index);
        end
    end
        %% Effective channel calculation
    TU_H_eff_2 = zeros(N_BS,Ns_total,N_sub);
    TU_W_HBF_fixed_2 = zeros(N_BS,Ns_total,N_sub);
    TU_result_CH_2 = zeros(Ns_total,N_sub);
    Noise_eff_2 = zeros(Ns_total,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            for nn_2 = 1:Ns_equal_alloc_2(u)  
            TU_W_HBF_fixed_2(:,sum(Ns_equal_alloc_2(1:1:u-1))+1:1:sum(Ns_equal_alloc_2(1:1:u)),sub_index) = TU_W_HBFu_fixed_2(:,1:Ns_equal_alloc_2(u),u,sub_index);
            TU_H_eff_2(:,sum(Ns_equal_alloc_2(1:1:u-1))+1:1:sum(Ns_equal_alloc_2(1:1:u)),sub_index) = Sel_H_user_2(:,:,sub_index,u)*TU_F_HBFu_fixed_2(:,1:Ns_equal_alloc_2(u),u,sub_index);
            
            Noise_eff_2(sum(Ns_equal_alloc_2(1:1:u-1))+nn_2,sub_index) = N0*norm(TU_W_HBF_fixed_2(sum(Ns_equal_alloc_2(1:1:u-1))+nn_2,sum(Ns_equal_alloc_2(1:1:u-1))+nn_2,sub_index)')^2;
            Noise_eff_2;
            TU_result_CH_2(sum(Ns_equal_alloc_2(1:1:u-1))+nn_2,sub_index) = norm(TU_W_HBF_fixed_2(sum(Ns_equal_alloc_2(1:1:u-1))+nn_2,sum(Ns_equal_alloc_2(1:1:u-1))+nn_2,sub_index)'*TU_H_eff_2(sum(Ns_equal_alloc_2(1:1:u-1))+nn_2,sum(Ns_equal_alloc_2(1:1:u-1))+nn_2,sub_index))^2;  %% Heff
            end
        end  
    end
    TU_result_CH_2;
         %% Step3. MODCOD index loading 
        Unit_bits = data_rate_change(1)*ones(Ns_total,N_sub);
        SINR = switching_level_change(1)*ones(Ns_total,N_sub); 
        MODCOD_Index_table3 = zeros(Ns_total,N_sub);                                                         % Set up MODCOD index table
        PT_unit_table3 = (SINR./Unit_bits.*Noise_eff_2)./TU_result_CH_2;   
        for nu = 1:U 
            st1 = sum(Ns_equal_alloc_2(1:(nu-1)))+1;
            ed1 = sum(Ns_equal_alloc_2(1:nu)); 
            bits = N_bits;                                                                             % Available bits                           
            while bits > 0
                [x1,y1] = find(PT_unit_table3(st1:ed1,:) == min(min(PT_unit_table3(st1:ed1,:))));                          % Find least power increase 
                MODCOD_Index_table3(st1+x1(1)-1,y1(1)) = MODCOD_Index_table3(st1+x1(1)-1,y1(1)) + 1;                             % Update MODCOD index table
                PT_unit_table3(st1+x1(1)-1,y1(1)) = (switching_level_change(MODCOD_Index_table3(st1+x1(1)-1,y1(1))+1)/data_rate_change(MODCOD_Index_table3(st1+x1(1)-1,y1(1))+1)*Noise_eff_2(st1+x1(1)-1,y1(1)))./TU_result_CH_2(st1+x1(1)-1,y1(1));
                bits = bits - data_rate_change(MODCOD_Index_table3(st1+x1(1)-1,y1(1)));                                    % Update transmitted power efficiency table                                                                                                                                                   % Update available bits
            end
        end 
        MODCOD_Index_table3;
        bits_table1 = data_rate_st(MODCOD_Index_table3+1);
        %% Calculate transmitted power
        for k = 1:N_sub
            for nf3 = 1:Ns_total
                PT_final_table3(nf3,k) = (switching_level_st(MODCOD_Index_table3(nf3,k)+1)*Noise_eff_2(nf3,k))/TU_result_CH_2(nf3,k);  
            end
        end 
        PT_final_table3;
        %% TUMD (K/4)
    % analog precoder at MS via Tensor unfolding
    long_unfold_Hu_4 = zeros((N_sub/4)*N_BS,N_per_MS,U); 
    long_unfold_F_RFu_4 = zeros(N_per_MS,Nrf_per_MS,U);
    for u=1:1:U
        for sub_index3=2:4:14
            long_unfold_Hu_4((sub_index3-1)*N_BS+1:1:sub_index3*N_BS,:,u) = Sel_H_user_2(:,:,sub_index3,u);
        end
        [long_unfold_V_4,long_unfold_D_4] = eig(long_unfold_Hu_4(:,:,u)'*long_unfold_Hu_4(:,:,u));  
        [~,longunf_ind_4] = sort(diag(long_unfold_D_4),'descend');
        long_unfold_F_RFu_4(:,:,u) = extract_phase_v2(long_unfold_V_4(:,longunf_ind_4(1:1:Nrf_per_MS)),N_per_MS);
    end    
    
    % analog combiner at BS via Tensor unfolding and RF chain assignment
    horiz_unfold_Hu_4 = zeros(N_BS,Nrf_per_MS*(N_sub/4),U);
    horiz_unfold_W_RF_4 = [];
    eigval_buff_4 = zeros(N_BS,U);
    eigvec_buff_4 = zeros(N_BS,N_BS,U);
    
    for u=1:1:U
        for sub_index3=2:4:14
            horiz_unfold_Hu_4(:,(sub_index3-1)*Nrf_per_MS+1:1:sub_index3*Nrf_per_MS,u) = Sel_H_user_2(:,:,sub_index3,u)*long_unfold_F_RFu_4(:,:,u);
        end
        [horiz_unfold_V_4,horiz_unfold_D_4] = eig(horiz_unfold_Hu_4(:,:,u)*horiz_unfold_Hu_4(:,:,u)');
        [horizunfout_4,horizunfind_4] = sort(diag(horiz_unfold_D_4),'descend');
        
        eigval_buff_4(:,u) = real(horizunfout_4);
        eigvec_buff_4(:,:,u) = horiz_unfold_V_4(:,horizunfind_3);  
    end   
 
    for u=1:1:U
        horiz_unfold_W_RF_4 = [horiz_unfold_W_RF_4,extract_phase_v2(eigvec_buff_4(:,1:1:Ns_TU_alloc_4(u),u),N_BS)];
    end
    
    % BD(ignore the noise)+SVD beamforming based on effective channel
    [TU_F_BBu_4,TU_W_BBu_4]=UPBB_general_BDSVD(long_unfold_F_RFu_4,horiz_unfold_W_RF_4,Sel_H_user_2,U,Ns_TU_alloc_4,N_sub,Nrf_per_MS,Nrf_BS);
    
    % Compute hybrid beamforming matrix
    TU_F_HBFu_4 = zeros(N_per_MS,max(Ns_TU_alloc_4),U,N_sub);
    TU_W_HBFu_4 = zeros(N_BS,max(Ns_TU_alloc_4),U,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            TU_F_HBFu_4(:,1:Ns_TU_alloc_4(u),u,sub_index) = long_unfold_F_RFu_4(:,:,u)*TU_F_BBu_4(:,1:Ns_TU_alloc_4(u),u,sub_index);
            
            % normalize transmit HBF Weight vector to one
            for ns = 1:1:Ns_TU_alloc_4(u)
                TU_F_HBFu_4(:,ns,u,sub_index) = TU_F_HBFu_4(:,ns,u,sub_index)/norm(TU_F_HBFu_4(:,ns,u,sub_index)); 
            end
            
            TU_W_HBFu_4(:,1:Ns_TU_alloc_4(u),u,sub_index) = horiz_unfold_W_RF_4*TU_W_BBu_4(:,1:Ns_TU_alloc_4(u),u,sub_index);
        end
    end
    %% Effective channel calculation
    TU_H_eff_4 = zeros(N_BS,Ns_total,N_sub);
    TU_W_HBF_fixed_4 = zeros(N_BS,Ns_total,N_sub);
    TU_result_CH_4 = zeros(Ns_total,N_sub);
    Noise_eff_4 = zeros(Ns_total,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            for nn_4 = 1:Ns_TU_alloc_4(u)
            TU_W_HBF_fixed_4(:,sum(Ns_TU_alloc_4(1:1:u-1))+1:1:sum(Ns_TU_alloc_4(1:1:u)),sub_index) = TU_W_HBFu_4(:,1:Ns_TU_alloc_4(u),u,sub_index);
            TU_H_eff_4(:,sum(Ns_TU_alloc_4(1:1:u-1))+1:1:sum(Ns_TU_alloc_4(1:1:u)),sub_index) = Sel_H_user_2(:,:,sub_index,u)*TU_F_HBFu_4(:,1:Ns_TU_alloc_4(u),u,sub_index);
            Noise_eff_4(sum(Ns_TU_alloc_4(1:1:u-1))+nn_4,sub_index) = N0*norm(TU_W_HBF_fixed_4(sum(Ns_TU_alloc_4(1:1:u-1))+nn_4,sum(Ns_TU_alloc_4(1:1:u-1))+nn_4,sub_index)')^2;
            Noise_eff_4;
            TU_result_CH_4(sum(Ns_TU_alloc_4(1:1:u-1))+nn_4,sub_index) = norm(TU_W_HBF_fixed_4(sum(Ns_TU_alloc_4(1:1:u-1))+nn_4,sum(Ns_TU_alloc_4(1:1:u-1))+nn_4,sub_index)'*TU_H_eff_4(sum(Ns_TU_alloc_4(1:1:u-1))+nn_4,sum(Ns_TU_alloc_4(1:1:u-1))+nn_4,sub_index))^2;  %% Heff
            end  
        end
    end
    TU_result_CH_4;
         %% Step3. MODCOD index loading 
        Unit_bits = data_rate_change(1)*ones(Ns_total,N_sub);
        SINR = switching_level_change(1)*ones(Ns_total,N_sub); 
        MODCOD_Index_table1_4 = zeros(Ns_total,N_sub);                                                         % Set up MODCOD index table
        PT_unit_table1_4 = (SINR./Unit_bits.*Noise_eff_4)./TU_result_CH_4;   
        for nu = 1:U 
            st = sum(Ns_TU_alloc_4(1:(nu-1)))+1;
            ed = sum(Ns_TU_alloc_4(1:nu)); 
            bits = N_bits;                                                                             % Available bits                           
            while bits > 0
                [x,y] = find(PT_unit_table1_4(st:ed,:) == min(min(PT_unit_table1_4(st:ed,:))));                          % Find least power increase 
                MODCOD_Index_table1_4(st+x(1)-1,y(1)) = MODCOD_Index_table1_4(st+x(1)-1,y(1)) + 1;                             % Update MODCOD index table
                PT_unit_table1_4(st+x(1)-1,y(1)) = (switching_level_change(MODCOD_Index_table1_4(st+x(1)-1,y(1))+1)/data_rate_change(MODCOD_Index_table1_4(st+x(1)-1,y(1))+1)*Noise_eff_4(st+x(1)-1,y(1)))./TU_result_CH_4(st+x(1)-1,y(1));
                bits = bits - data_rate_change(MODCOD_Index_table1_4(st+x(1)-1,y(1)));                                    % Update transmitted power efficiency table                                                                                                                                                   % Update available bits
            end
        end 
        MODCOD_Index_table1_4;
        bits_table1_4 = data_rate_st(MODCOD_Index_table1_4+1);
        %% Calculate transmitted power
        for k = 1:N_sub
            for nf = 1:Ns_total
                PT_final_table1_4(nf,k) = (switching_level_st(MODCOD_Index_table1_4(nf,k)+1)*Noise_eff_4(nf,k))/TU_result_CH_4(nf,k);  
            end
        end 
        PT_final_table1_4;
        %% Fully-digital beamforming
        FD_Hu_2 = zeros(N_BS,N_per_MS*(U),N_sub); 
        FD_U_tilte_2 = zeros(N_BS,N_BS,U,N_sub);
        rank_int_ch_2 = zeros(U,N_sub);
    for sub_index=1:1:N_sub
        for u=1:1:U
            FD_Hu_2(:,(u-1)*N_per_MS+1:1:u*N_per_MS,sub_index) = Sel_H_user_2(:,:,sub_index,u);
        end
        for u=1:1:U
            FD_Hu_tilte_2 = FD_Hu_2(:,:,sub_index);
            FD_Hu_tilte_2(:,(u-1)*N_per_MS+1:1:u*N_per_MS) = [];
            
            [FD_U_tilte_2(:,:,u,sub_index),~,~] = svd(FD_Hu_tilte_2);
            rank_int_ch_2(u,sub_index) = rank(FD_Hu_tilte_2);
        end     
    end

        FD_W_u_2 = zeros(N_BS,max(Ns_FD_alloc_2),U,N_sub);
        FD_F_u_2 = zeros(N_per_MS,max(Ns_FD_alloc_2),U,N_sub);
    for sub_index=1:1:N_sub
        for u=1:1:U
            [FD_UU_tilte_2,~,FD_VV_tilte_2] = svd(FD_U_tilte_2(:,rank_int_ch_2(u,sub_index)+1:end,u,sub_index)'*Sel_H_user_2(:,:,sub_index,u));
            FD_F_u_2(:,1:Ns_FD_alloc_2(u),u,sub_index) = FD_VV_tilte_2(:,1:Ns_FD_alloc_2(u));
            FD_W_u_2(:,1:Ns_FD_alloc_2(u),u,sub_index) = FD_U_tilte_2(:,rank_int_ch_2(u,sub_index)+1:end,u,sub_index)*FD_UU_tilte_2(:,1:Ns_FD_alloc_2(u));    
        end     
    end
        %% Effective channel calculation
    FD_H_eff_2 = zeros(N_BS,Ns_total,N_sub);
    FD_W_HBF_fixed_2 = zeros(N_BS,Ns_total,N_sub);
    FD_result_CH_2 = zeros(Ns_total,N_sub);
    FD_Noise_eff_2 = zeros(Ns_total,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            for fn_2 = 1:Ns_FD_alloc_2(u)
            FD_W_HBF_fixed_2(:,sum(Ns_FD_alloc_2(1:1:u-1))+1:1:sum(Ns_FD_alloc_2(1:1:u)),sub_index) = FD_W_u_2(:,1:Ns_FD_alloc_2(u),u,sub_index);
            FD_H_eff_2(:,sum(Ns_FD_alloc_2(1:1:u-1))+1:1:sum(Ns_FD_alloc_2(1:1:u)),sub_index) = Sel_H_user_2(:,:,sub_index,u)*FD_F_u_2(:,1:Ns_FD_alloc_2(u),u,sub_index);
            
            FD_Noise_eff_2(sum(Ns_FD_alloc_2(1:1:u-1))+fn_2,sub_index) = N0*norm(FD_W_HBF_fixed_2(sum(Ns_FD_alloc_2(1:1:u-1))+fn_2,sum(Ns_FD_alloc_2(1:1:u-1))+fn_2,sub_index)')^2;
            FD_result_CH_2(sum(Ns_FD_alloc_2(1:1:u-1))+fn_2,sub_index) = norm(FD_W_HBF_fixed_2(sum(Ns_FD_alloc_2(1:1:u-1))+fn_2,sum(Ns_FD_alloc_2(1:1:u-1))+fn_2,sub_index)'*FD_H_eff_2(sum(Ns_FD_alloc_2(1:1:u-1))+fn_2,sum(Ns_FD_alloc_2(1:1:u-1))+fn_2,sub_index))^2;  %% Heff
            end
        end
    end
    FD_result_CH_2;
         %% Step3. MODCOD index loading 
        Unit_bits = data_rate_change(1)*ones(Ns_total,N_sub);
        SINR = switching_level_change(1)*ones(Ns_total,N_sub); 
        MODCOD_Index_table4 = zeros(Ns_total,N_sub);                                                         % Set up MODCOD index table
        PT_unit_table4 = (SINR./Unit_bits.*FD_Noise_eff_2)./FD_result_CH_2;                 % Set up transmitted power efficiency table    
        for nu = 1:U
            st1 = sum(Ns_FD_alloc_2(1:1:nu-1))+1;
            ed1 = sum(Ns_FD_alloc_2(1:1:nu)); 
            bits = N_bits;                                                                             % Available bits                           
            while bits > 0
                [x1,y1] = find(PT_unit_table4(st1:ed1,:) == min(min(PT_unit_table4(st1:ed1,:))));                          % Find least power increase 
       
                MODCOD_Index_table4(st1+x1(1)-1,y1(1)) = MODCOD_Index_table4(st1+x1(1)-1,y1(1)) + 1;                             % Update MODCOD index table
                PT_unit_table4(st1+x1(1)-1,y1(1)) = (switching_level_change(MODCOD_Index_table4(st1+x1(1)-1,y1(1))+1)/data_rate_change(MODCOD_Index_table4(st1+x1(1)-1,y1(1))+1)*FD_Noise_eff_2(st1+x1(1)-1,y1(1)))./FD_result_CH_2(st1+x1(1)-1,y1(1));
                bits = bits - data_rate_change(MODCOD_Index_table4(st1+x1(1)-1,y1(1)));                                    % Update transmitted power efficiency table                                                                                                                                                   % Update available bits
            end
        end 
        MODCOD_Index_table4;
        bits_table2 = data_rate_st(MODCOD_Index_table4+1);
        %% Calculate transmitted power
        for k = 1:N_sub
            for nf4 = 1:Ns_total
                PT_final_table4(nf4,k) = (switching_level_st(MODCOD_Index_table4(nf4,k)+1)*FD_Noise_eff_2(nf4,k))/FD_result_CH_2(nf4,k);
            end
        end 
        PT_final_table4;
        %% Benchmark Scheme--ACMD
        % analog precoder at MS and analog combiner at BS via channel average
        chavg_F_RFu_2 = zeros(N_per_MS,Nrf_per_MS,U);
        chavg_W_RF_2 = [];
        chavg_sinval_buff_2 = zeros(N_per_MS,U);
        chavg_L_sinvec_buff_2 = zeros(N_BS,N_per_MS,U);
        % analog precoder
        for u=1:1:U
            Q_2 = 0;
            for sub_index=1:1:N_sub
                Q_2 = Q_2+Sel_H_user_2(:,:,sub_index,u)/N_sub;
            end
            [chavg_U_2,~,chavg_V_2] = svd(Q_2);
            chavg_F_RFu_2(:,:,u) = extract_phase_v2(chavg_V_2(:,1:Nrf_per_MS),N_per_MS);
            chavg_sinval_buff_2(:,u) = svd(Q_2);
            chavg_L_sinvec_buff_2(:,:,u) = chavg_U_2(:,1:N_per_MS);
        end  

        for u=1:1:U
            chavg_W_RF_2 = [chavg_W_RF_2,extract_phase_v2(chavg_L_sinvec_buff_2(:,1:1:Ns_CHAVG_alloc_2(u),u),N_BS)];
        end

        % BD(ignore the noise)+SVD beamforming based on effective channel
        [chavg_F_BBu_2,chavg_W_BBu_2]=UPBB_general_BDSVD(chavg_F_RFu_2,chavg_W_RF_2,Sel_H_user_2,U,Ns_CHAVG_alloc_2,N_sub,Nrf_per_MS,Nrf_BS);

        % Compute hybrid beamforming matrix
        chavg_F_HBFu_2 = zeros(N_per_MS,max(Ns_CHAVG_alloc_2),U,N_sub);
        chavg_W_HBFu_2 = zeros(N_BS,max(Ns_CHAVG_alloc_2),U,N_sub);
        for sub_index=1:1:N_sub       
            for u=1:1:U
                chavg_F_HBFu_2(:,1:Ns_CHAVG_alloc_2(u),u,sub_index) = chavg_F_RFu_2(:,:,u)*chavg_F_BBu_2(:,1:Ns_CHAVG_alloc_2(u),u,sub_index);

                % normalize transmit HBF Weight vector to one
                for ns = 1:1:Ns_CHAVG_alloc_2(u)
                    chavg_F_HBFu_2(:,ns,u,sub_index) = chavg_F_HBFu_2(:,ns,u,sub_index)/norm(chavg_F_HBFu_2(:,ns,u,sub_index)); 
                end

                chavg_W_HBFu_2(:,1:Ns_CHAVG_alloc_2(u),u,sub_index) = chavg_W_RF_2*chavg_W_BBu_2(:,1:Ns_CHAVG_alloc_2(u),u,sub_index);
            end
        end
    %% Effective channel calculation
    chavg_H_eff_2 = zeros(N_BS,Ns_total,N_sub);
    chavg_W_HBF_fixed_2 = zeros(N_BS,Ns_total,N_sub);
    chavg_result_CH_2 = zeros(Ns_total,N_sub);
    chavg_Noise_eff_2 = zeros(Ns_total,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            for cn_2 = 1:Ns_CHAVG_alloc_2(u)
            chavg_W_HBF_fixed_2(:,sum(Ns_CHAVG_alloc_2(1:1:u-1))+1:1:sum(Ns_CHAVG_alloc_2(1:1:u)),sub_index) = chavg_W_HBFu_2(:,1:Ns_CHAVG_alloc_2(u),u,sub_index);
            chavg_H_eff_2(:,sum(Ns_CHAVG_alloc_2(1:1:u-1))+1:1:sum(Ns_CHAVG_alloc_2(1:1:u)),sub_index) = Sel_H_user_2(:,:,sub_index,u)*chavg_F_HBFu_2(:,1:Ns_CHAVG_alloc_2(u),u,sub_index);
            
            chavg_Noise_eff_2(sum(Ns_CHAVG_alloc_2(1:1:u-1))+cn_2,sub_index) = N0*norm(chavg_W_HBF_fixed_2(sum(Ns_CHAVG_alloc_2(1:1:u-1))+cn_2,sum(Ns_CHAVG_alloc_2(1:1:u-1))+cn_2,sub_index)')^2;
            chavg_result_CH_2(sum(Ns_CHAVG_alloc_2(1:1:u-1))+cn_2,sub_index) = norm(chavg_W_HBF_fixed_2(sum(Ns_CHAVG_alloc_2(1:1:u-1))+cn_2,sum(Ns_CHAVG_alloc_2(1:1:u-1))+cn_2,sub_index)'*chavg_H_eff_2(sum(Ns_CHAVG_alloc_2(1:1:u-1))+cn_2,sum(Ns_CHAVG_alloc_2(1:1:u-1))+cn_2,sub_index))^2;  %% Heff
            end
        end
    end
    chavg_result_CH_2;
         %% Step3. MODCOD index loading 
        Unit_bits = data_rate_change(1)*ones(Ns_total,N_sub);
        SINR = switching_level_change(1)*ones(Ns_total,N_sub); 
        MODCOD_Index_table6 = zeros(Ns_total,N_sub);                                                         % Set up MODCOD index table
        PT_unit_table6 = (SINR./Unit_bits.*chavg_Noise_eff_2)./chavg_result_CH_2;                 % Set up transmitted power efficiency table    
        for nu = 1:U
            st1 = sum(Ns_CHAVG_alloc_2(1:1:nu-1))+1;
            ed1 = sum(Ns_CHAVG_alloc_2(1:1:nu)); 
            bits = N_bits;                                                                             % Available bits                           
            while bits > 0
                [x1,y1] = find(PT_unit_table6(st1:ed1,:) == min(min(PT_unit_table6(st1:ed1,:))));                          % Find least power increase 
                MODCOD_Index_table6(st1+x1(1)-1,y1(1)) = MODCOD_Index_table6(st1+x1(1)-1,y1(1)) + 1;                             % Update MODCOD index table
                PT_unit_table6(st1+x1(1)-1,y1(1)) = (switching_level_change(MODCOD_Index_table6(st1+x1(1)-1,y1(1))+1)/data_rate_change(MODCOD_Index_table6(st1+x1(1)-1,y1(1))+1)*chavg_Noise_eff_2(st1+x1(1)-1,y1(1)))./chavg_result_CH_2(st1+x1(1)-1,y1(1));
                bits = bits - data_rate_change(MODCOD_Index_table6(st1+x1(1)-1,y1(1)));                                    % Update transmitted power efficiency table                                                                                                                                                   % Update available bits
            end
        end 
        MODCOD_Index_table6;
        bits_table6 = data_rate_st(MODCOD_Index_table6+1);
        %% Calculate transmitted power
        for k = 1:N_sub
            for nf1 = 1:Ns_total
                PT_final_table6(nf1,k) = (switching_level_st(MODCOD_Index_table6(nf1,k)+1)*chavg_Noise_eff_2(nf1,k))/chavg_result_CH_2(nf1,k);
            end
        end 
        PT_final_table6;

        PT_final1(ch,r) = sum(sum(PT_final_table1));
        PT_final1_3(ch,r) = sum(sum(PT_final_table1_3));
        PT_final1_4(ch,r) = sum(sum(PT_final_table1_4));
        PT_final2(ch,r) = sum(sum(PT_final_table2));
        PT_final3(ch,r) = sum(sum(PT_final_table3));
        PT_final4(ch,r) = sum(sum(PT_final_table4));
        PT_final5(ch,r) = sum(sum(PT_final_table5));
        PT_final6(ch,r) = sum(sum(PT_final_table6));
    end % end of U
end % end of channel realization

if ch_reali == 1
    PT_final_mean1 = PT_final1;
    PT_final_mean1_3 = PT_final1_3;
    PT_final_mean1_4 = PT_final1_4;
    PT_final_mean2 = PT_final2;
    PT_final_mean3 = PT_final3;
    PT_final_mean4 = PT_final4;
    PT_final_mean5 = PT_final5;
    PT_final_mean6 = PT_final6;
else
    PT_final_mean1 = sum(PT_final1)/ch_reali;
    PT_final_mean1_3 = sum(PT_final1_3)/ch_reali;
    PT_final_mean1_4 = sum(PT_final1_4)/ch_reali;
    PT_final_mean2 = sum(PT_final2)/ch_reali;
    PT_final_mean3 = sum(PT_final3)/ch_reali;
    PT_final_mean4 = sum(PT_final4)/ch_reali;
    PT_final_mean5 = sum(PT_final5)/ch_reali;
    PT_final_mean6 = sum(PT_final6)/ch_reali;
end
norm_power1 = 10*log10(PT_final_mean1/PT_final_mean6(1));
norm_power1_3 = 10*log10(PT_final_mean1_3/PT_final_mean6(1));
norm_power1_4 = 10*log10(PT_final_mean1_4/PT_final_mean6(1));
norm_power2 = 10*log10(PT_final_mean2/PT_final_mean6(1));
norm_power3 = 10*log10(PT_final_mean3/PT_final_mean6(1));
norm_power4 = 10*log10(PT_final_mean4/PT_final_mean6(1));
norm_power5 = 10*log10(PT_final_mean5/PT_final_mean6(1));
norm_power6 = 10*log10(PT_final_mean6/PT_final_mean6(1));

figure
plot(Total_U,norm_power1,'bo-','LineWidth',1.2,'DisplayName','Proposed (US)')
hold on
plot(Total_U,norm_power1_3,'mo-.','LineWidth',1.2,'DisplayName','Proposed (Reduce Complexity,K/4)(US)')
plot(Total_U,norm_power2,'-*','color',[0.1500 0.8000 0.3000],'LineWidth',1.2,'DisplayName','Fully-Digital (US)')
plot(Total_U,norm_power5,'-^','color',[0.8500 0.3250 0.0980],'LineWidth',1.2,'DisplayName','Modified\_ACMD [35] (US)')

plot(Total_U,norm_power3,'rs-','LineWidth',1.2,'DisplayName','Proposed (RUS)')
plot(Total_U,norm_power1_4,'bs-.','LineWidth',1.2,'DisplayName','Proposed (Reduce Complexity,K/4)(RUS)')
plot(Total_U,norm_power4,'-.*','color',[0.1500 0.8000 0.3000],'LineWidth',1.2,'DisplayName','Fully-Digital (RUS)')
plot(Total_U,norm_power6,'-.^','color',[0.8500 0.3250 0.0980],'LineWidth',1.2,'DisplayName','Modified\_ACMD [35] (RUS)')

legend('Location','best')
xlabel('Number of MSs (U)')
ylabel('Relative Transmitted Power (dB)')
% axis([-inf inf -inf inf])
grid on
hold off
set(gca,'XTick',Total_U)
toc