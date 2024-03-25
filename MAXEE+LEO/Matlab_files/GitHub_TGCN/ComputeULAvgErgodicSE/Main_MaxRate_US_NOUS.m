% Upink multi-ant
tic;clc;clear;close all;
%% parameter seting
% BS parameter
N_BS=128;
Nrf_BS=8;                                                                  % number of RF chanis at BS(This parameter determine the most number of receive data stream)

% Mobile parameter
% ?™è£¡??„Total_U?‚ºcell?…§ç¸½å…±??„ä½¿?”¨?…æ•¸ï¼Œä½¿?”¨?…é¸??‡æ?Ÿåˆ¶è¦å?ä¸­?¸U?‹ä½¿?”¨?…å‡ºä¾?
Total_U = 20;                                                              
U = 4; 

Ns_per_MS=2;                                                               % number of data streams per user
Nrf_per_MS=3;                                                              % number of RF chanis per user(This parameter determine the most number of transmitt data stream per user) 
N_per_MS=16;
if Nrf_BS<=Nrf_per_MS*U
    Ns_total = Nrf_BS;
else
    Ns_total = Nrf_per_MS*U;
end
Ns_TU_alloc = Ns_per_MS*ones(1,U);
Ns_FD_alloc = Ns_per_MS*ones(1,U);
Ns_CHAVG_alloc = Ns_per_MS*ones(1,U);
Ns_TU_alloc_2 = Ns_per_MS*ones(1,U);
Ns_FD_alloc_2 = Ns_per_MS*ones(1,U);
Ns_CHAVG_alloc_2 = Ns_per_MS*ones(1,U);

h_BS=10;                                                                   % BS height[meter]
R=100;                                                                     % distance of BS and user[meter]
d0=10;                                                                     % minimum  distance between BS and UE(meter)
fc=28;                                                                     % carrier frequency(28 GHz)
N_sub=16;                                                                  % number of subcarriers

System_BW=(20*10^6);                                                       % System Bandwidth(50 MHz)
sub_BW=System_BW/N_sub;                                                    % Subcarrier frequency spacing
noise_variance_dBm=-159+10*log10(sub_BW);                                  % noise density=-159 dBm/Hz
N0=10^(noise_variance_dBm/10)*10^-3;                                       % noise power[W] at each subcarrier

% Ref:5G;NR; User Equipment (UE) radio transmission and reception;Part 2:
% Range 2 Standalone (3GPP TS 38.101-2 version 15.3.0 Release 15) CH 6.
UE_Max_Tx_dBm=30;                                                          % UE maximum output power limits for power class 3(dBm)
P_Max=(10^(UE_Max_Tx_dBm/10)*10^-3);                                       % UE maximum output power limits for power class 3(W) on entire bandwidth
UE_TX_POW_dBm=5:5:20;
ch_reali=10^8;
P_total_user=(10.^(UE_TX_POW_dBm/10)/1000);

                                                                                                                       
ang_spread=10*pi/180;
Nray=10;
%% main program
for ch=1:1:ch_reali
       %% UE are randomly dropped in the cell over distances ranging between 10 m and the cell radius
    % ?™é?Šè?æ?ŠTotal_U?‹ä½¿?”¨?…å?„è‡ª?œ¨cell??„ä?ç½®?”¢??Ÿå‡ºä¾†ï?ŒU?”¹??Total_U
    user_location_phase=pi*(rand(1,Total_U)*2-1);%-pi~pi
    user_location_radius=rand(1,Total_U)*(R-d0)+d0; 
    user_location=user_location_radius.*exp(j*user_location_phase);
    x_location=real(user_location);
    y_location=imag(user_location);
       %% clustered channel model
    H_user = zeros(N_BS,N_per_MS,N_sub,Total_U);
    % ?™å?‹forè¿´å?ˆä?æ¨?è¦ç?™è?—ï?Œå? ç‚ºè¦å?ˆæ?ŠTotal_U?‹ä½¿?”¨?…å?„è‡ª??„é?šé?“å?ˆç”¢??Ÿå‡ºä¾?
    for u=1:1:Total_U
        d_2D(u)=norm([x_location(u) y_location(u)],2);                     % unit is meter    
        [PL(u),path_loss_dB(u)]=PathLoss(fc,h_BS,d_2D(u));                            % path loss model
        [CH,Ar,At,Sort_beta,~] =Satellite_channel(N_per_MS,N_BS,N_sub,1,Nray,ang_spread,path_loss_dB(u),fc);
        H_user(:,:,:,u)=CH(:,:,:,1)/sqrt(PL(u));
    end
       %% TUMD
    % analog precoder at MS via Tensor unfolding
    long_unfold_Hu = zeros((N_sub)*N_BS,N_per_MS,U); 
    long_unfold_F_RFu = zeros(N_per_MS,Nrf_per_MS,U);
    for u=1:1:U
        for sub_index=1:1:N_sub
            % H_user ??›æ?? Sel_H_user
            long_unfold_Hu((sub_index-1)*N_BS+1:1:sub_index*N_BS,:,u) = H_user(:,:,sub_index,u);
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
            % H_user ??›æ?? Sel_H_user
            horiz_unfold_Hu(:,(sub_index-1)*Nrf_per_MS+1:1:sub_index*Nrf_per_MS,u) = H_user(:,:,sub_index,u)*long_unfold_F_RFu(:,:,u);
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
    % H_user ??›æ?? Sel_H_user
    [TU_F_BBu,TU_W_BBu]=UPBB_general_BDSVD(long_unfold_F_RFu,horiz_unfold_W_RF,H_user,U,Ns_TU_alloc,N_sub,Nrf_per_MS,Nrf_BS);
    
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
       %% Scheme 1:  Proposed User Selection
           max_userNs=4;
           max_singular = zeros(N_sub,Total_U,1);
        for nnu=1:Total_U
           for nc = 1:N_sub
              [~,sl(:,:,nc),~] = svd(H_user(:,:,nc,nnu),'econ'); 
              sl_(:,:,nc) = diag(sl(:,:,nc));               
           end
           for chh = 1:1
               max_singular(:,nnu,chh) = reshape(sl_(chh,1,:),1,N_sub);
               [sorted_gain(:,:,chh),IX(:,:,chh)] = sort(max_singular(:,:,chh),2,'descend');
           end          
        end

        selected_subcarriers = IX(:,1:U);

% find maximum channel        
if max_userNs == 1
    for nc = 1:N_sub
        for nuu = 1:Total_U
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
        for nu1 = 1:Total_U
            if selected_subcarriers(nsc,nu1)~=0
                channel_max(nsc,i) = selected_subcarriers(nsc,nu1);
                i=i+1;
            end
            if i == max_userNs+1       % Ns=2 - i=2
                break;
            end
        end
    end 
end
selected_subcarriers = sort(selected_subcarriers,2,'descend');    
Sel_User_ind_2 = selected_subcarriers(1,1:U);
       %% ?™é?Šæ˜¯?‚ºäº†æ?Šè¢«?¸?ˆ°??„ä½¿?”¨?…é?šé?“æ?‘å‡ºä¾†å?˜æ”¾
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
    
    % analog combiner at BS via Tensor unfolding and RF chain assignment
    horiz_unfold_Hu_2 = zeros(N_BS,Nrf_per_MS*(N_sub),U);
    horiz_unfold_W_RF_2 = [];
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
 
    for u=1:1:U
        horiz_unfold_W_RF_2 = [horiz_unfold_W_RF_2,extract_phase_v2(eigvec_buff_2(:,1:1:Ns_TU_alloc_2(u),u),N_BS)];
    end
    
    % BD(ignore the noise)+SVD beamforming based on effective channel
    [TU_F_BBu_2,TU_W_BBu_2]=UPBB_general_BDSVD(long_unfold_F_RFu_2,horiz_unfold_W_RF_2,Sel_H_user_2,U,Ns_TU_alloc_2,N_sub,Nrf_per_MS,Nrf_BS);
    
    % Compute hybrid beamforming matrix
    TU_F_HBFu_2 = zeros(N_per_MS,max(Ns_TU_alloc_2),U,N_sub);
    TU_W_HBFu_2 = zeros(N_BS,max(Ns_TU_alloc_2),U,N_sub);
    for sub_index=1:1:N_sub       
        for u=1:1:U
            TU_F_HBFu_2(:,1:Ns_TU_alloc_2(u),u,sub_index) = long_unfold_F_RFu_2(:,:,u)*TU_F_BBu_2(:,1:Ns_TU_alloc_2(u),u,sub_index);
            
            % normalize transmit HBF Weight vector to one
            for ns = 1:1:Ns_TU_alloc_2(u)
                TU_F_HBFu_2(:,ns,u,sub_index) = TU_F_HBFu_2(:,ns,u,sub_index)/norm(TU_F_HBFu_2(:,ns,u,sub_index)); 
            end
            
            TU_W_HBFu_2(:,1:Ns_TU_alloc_2(u),u,sub_index) = horiz_unfold_W_RF_2*TU_W_BBu_2(:,1:Ns_TU_alloc_2(u),u,sub_index);
        end
    end
       %%
    for r=1:1:length(UE_TX_POW_dBm)
        [ch r]
             %% -------------------------------Water filling power allocation ------------------------- 
        % Scheme 1:      
        % H_user ??›æ?? Sel_H_user
        [WFPow_TU,CNR_eff_TU] = Multi_Ns_WFPA(P_total_user(r),TU_F_HBFu,TU_W_HBFu,H_user,U,Ns_TU_alloc,N_sub,N0);
        % Scheme 2:      
        [WFPow_TU_2,CNR_eff_TU_2] = Multi_Ns_WFPA(P_total_user(r),TU_F_HBFu_2,TU_W_HBFu_2,Sel_H_user_2,U,Ns_TU_alloc_2,N_sub,N0);
        
             %%  -------------------------------performance evaluation(multi user)--------------------    
        % Scheme 1: 
        % H_user ??›æ?? Sel_H_user
        [R_TU_WFPow,R_iku_TU_WFPow,S_TU_WFPow,inter_stream_I_TU_WFPow,inter_user_I_TU_WFPow,I_TU_WFPow,N_TU_WFPow] = Multi_Ns_UP_rate(WFPow_TU,TU_F_HBFu,TU_W_HBFu,H_user,N0,U,Ns_TU_alloc,N_sub);        
        % Scheme 2: 
        [R_TU_WFPow_2,R_iku_TU_WFPow_2,S_TU_WFPow_2,inter_stream_I_TU_WFPow_2,inter_user_I_TU_WFPow_2,I_TU_WFPow_2,N_TU_WFPow_2] = Multi_Ns_UP_rate(WFPow_TU_2,TU_F_HBFu_2,TU_W_HBFu_2,Sel_H_user_2,N0,U,Ns_TU_alloc_2,N_sub);         

        
        fixed_ch_R_TU_WFPow(ch,r) = R_TU_WFPow;
        fixed_ch_R_TU_WFPow_2(ch,r) = R_TU_WFPow_2;
    end % end           
end % end of channel realization 

avg_rate_TU_WFPow = mean(fixed_ch_R_TU_WFPow,1);
avg_rate_TU_WFPow_2 = mean(fixed_ch_R_TU_WFPow_2,1);

figure
plot(UE_TX_POW_dBm,avg_rate_TU_WFPow_2,'bo-','LineWidth',1.2,'DisplayName','Proposed (With US)')
hold on

plot(UE_TX_POW_dBm,avg_rate_TU_WFPow,'ro--','LineWidth',1.2,'DisplayName','Proposed (Without US)')


legend('Location','best')
xlabel('Transmit Power per MS (dBm)')
ylabel('Spectral Efficiency(bits/s/Hz)')
axis([-inf inf -inf inf])
grid on
hold off
toc