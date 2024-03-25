%This Matlab script is the main of the article:
%
% Andrea Pizzo, Daniel Verenzuela, Luca Sanguinetti and Emil Bj繹rnson, "Network Deployment for Maximal Energy Efficiency 
% in Uplink with Multislope Path Loss," IEEE Transactions on Green Communications and Networking, Submitted to.
%
%This is version 1.0 (Last edited: 10-April-2018)
%
%License: This code is licensed under the GPLv2 license. If you in any way
%use this code for research that results in publications, please cite our
%original article listed above.
%
%
%This script compute the lower bound on the average ergodic spectral
%efficiency as in Theorem 2 (Section 3) of the article cited above. 
%This is achieved by calling 3 functions: 
%functionExampleSetup, functionChannelEstimates and functionComputeSE_UL. The function "functionExampleSetup" 
%deploys the UEs and the BSs randomly wihin the coverage area as described
%in Section 2. The function "functionChannelEstimates" computes all the
%uplink channel estimate for all UEs in the entire network. The function
%"functionComputeSE_UL" uses the true channel and the channel estimates to
%compute the uplink average ergodic spectral efficiency.


%Empty workspace and close figures
close all;
clear;
clc;

% 定義 satellite_channel 函數的參數
Nt_user = 100;            % 每個用戶的發射天線數量
Nr_BS = 128;              % 基站的接收天線數量
no_of_subcarriers = 16;   % 子載波的數量
no_of_users = 6;          % 用戶的數量
Nray = 10;                % 射線或散射體的數量（同時到達）
angspread = pi/180*10;    % 角度展開（以弧度為單位）
path_loss_dB = 100;       % 路徑損耗（示例分貝值）
fc = 2e9;                 % 載波頻率（例如，2 GHz）

% 使用定義的參數調用 satellite_channel 函數
[H_user] = Satellite_channel(Nt_user,Nr_BS,no_of_subcarriers,no_of_users,Nray,angspread,path_loss_dB,fc);

%Number of UEs per BS
K = 10;

%Define the range of BS antennas
Mrange = 100;

%Extract maximum number of BS antennas
Mmax = max(Mrange);

%Define the average pilot reuse factor
fRange = 1:1:8; 

%Select the number of setups with random UE locations
nbrOfSetups = 1;

%Select the number of channel realizations per setup
nbrOfRealizations = 1;

% Generates BS locations over (1) a fixed grid or
% (2) as a Homogeneous Poisson Point Process H-PPP
deployment = 2;

%指平面上（二維平面）從基站（BS）到用戶設備（UE）的水平距離(KM)
d_2D = 1000;

%表基地站（BS）的高度(Meter)
h_BS = 10;

%在這段程式碼中，deployment 變數用來決定基地台的部署方式，有兩種選擇：
%固定網格（Fixed Grid）：
%當 deployment 變數的值為 1 時，表示基底站在一個固定的網格上進行部署。 這意味著基地台的位置被預先確定，並且通常以規則的網格狀方式分佈在所選區域內。
%均勻泊松點過程（Homogeneous Poisson Point Process，H-PPP）：

%當 deployment 變數的值為 2 時，表示基地台依照均勻泊松點程序進行部署。 在這種部署方式下，基地台的位置是隨機的，並且遵循泊松點過程的統計特性，即在給定區域內的任意位置的基地台數量服從泊松分佈。
%根據模擬的需要和實際情況，可以選擇合適的基地台部署方式。 固定網格適用於需要控制基地台位置並具有規律性的場景，而H-PPP適用於更真實和隨機的場景。


% Generates the large scale path-loss coefficient using (1) single slope or
% (2) multi slope model
pathloss = 1;

% (1) Correlated or (2) Uncorellated channel fading
c = 1;

% Density of BSs for the random deployment
lambdaRange = fliplr([1 2:2:10]);

%這段程式碼是用來定義隨機部署中基站（BS）的密度。在這個情境下，`lambdaRange` 是一個向量，它包含了不同的基站密度值，這些值會在仿真中使用。
%在這個特定的例子中，`fliplr([1 2:2:8])` 的意思是建立一個向量，該向量包含了從 1 開始的奇數數列 `[1 3 5 ...]`，然後使用 `fliplr` 函式將這個向量反轉，所以最後的 `lambdaRange` 就變成了 `[1 1 1]`。
%這表示在這個仿真中，基站密度是固定的且為 1，即每平方公里有一個基站。這樣的設置可能用於在密集型網絡下進行效能評估，其中每個基站的服務範圍被合理地重疊，以提供充足的覆蓋和容量。

% Check the lengths of the variables
disp(length(lambdaRange));
disp(length(fRange));
disp(length(Mrange));

%Load fixed Propagation parameters %這段是呼叫SetPropagationParameters這隻程式碼
SetPropagationParameters;
PathLoss(fc,h_BS,d_2D);
%%

%Prepare to save simulation results
sumSE_MR = zeros(length(lambdaRange),length(fRange),length(Mrange),nbrOfSetups);
sumSE_ZF = zeros(length(lambdaRange),length(fRange),length(Mrange),nbrOfSetups);
sumSE_RZF = zeros(length(lambdaRange),length(fRange),length(Mrange),nbrOfSetups);
sumSE_SMMSE = zeros(length(lambdaRange),length(fRange),length(Mrange),nbrOfSetups);
sumSE_MMMSE = zeros(length(lambdaRange),length(fRange),length(Mrange),nbrOfSetups);

%這段程式碼保存到Allen_simulation中
% Check the dimensions of the variables
disp(size(sumSE_MMMSE));
disp(size(sumSE_SMMSE));
disp(size(sumSE_RZF));
disp(size(sumSE_ZF));
disp(size(sumSE_MR));
%% Go through all lambda values
for ll = 1:length(lambdaRange)
    
    %Output simulation progress
    disp([num2str(ll) ' lambda values out of ' num2str(length(lambdaRange))]);
    
% 這段 MATLAB 使用了一個 for 循環來遍歷 lambdaRange 向量,並輸出當前的模擬進度。允許我為您詳細解釋:
%1. for ll = 1:length(lambdaRange) 定義了一個循還,變量 ll 從 1 迭代到 lambdaRange 的長度。
%2. 在每次迭代中,ll 對應 lambdaRange 中的一個基站密度值。
%3. disp() 用於在 MATLAB 命令視窗輸出信息。
%4. num2str() 函數可以把數字轉換為字符串。
%5. \[\] 用於字符串連接,實現類似字符串格式化的效果。
%6. 輸出的字符串代表當前是第 ll 個 lambda 值,總共有 length(lambdaRange) 個。
%整體意思就是,在循環迭代 lambdaRange 中每一個基站密度時,輸出當前的模擬進度。這樣可以方便追蹤長時間運行的模擬實驗的進展情況。
%這種 Simple 的進度輸出在長時間運行的實驗中是非常有用的,既可以查看當前狀態,也可以估計剩余時間。總體上很好地體現了 MATLAB 在工程應用上輸出信息和傳達狀態的便利性。
    % Go through all setups
    for n = 1:nbrOfSetups
        
        %Output simulation progress
        disp([num2str(n) ' setups out of ' num2str(nbrOfSetups)]);
        
        if deployment  == 1
            
            L = 16;
            
        elseif deployment == 2
            
            %Generate the number of BSs in the area per km^2
            L = round(poissrnd(lambdaRange(ll)*(squareLength/1000)^2));
            
            %If the number is zero, then make a new realization
            while L == 0
                L = round(poissrnd(lambdaRange(ll)*(squareLength/1000)^2));
            end
    
        end
%這段程式碼是一個迴圈，用於遍歷模擬中的不同設定。在每個設定中，它會檢查部署的類型是否為 1 或 2。根據部署的類型，它會設置參數 L 的值。
%如果部署類型是 1，則 L 被設置為 16。這意味著在這種情況下，系統中有 16 個基站（或細胞）。
%如果部署類型是 2，則它使用 Poisson 分佈生成一個隨機的 BS 數量 L，這個數量是根據某個範圍內的平均 BS 數密度和方形區域的大小計算得出的。生成的 L 值可能不是整數，因此使用 round 函數進行四捨五入。同時，它確保 L 不會為零，如果生成的值為零，則重新生成。
%總之，這段程式碼確保在每個模擬設定中，基站（或細胞）的數量 L 都會根據部署類型而確定，這是模擬中的一個重要參數。


        squareLength = 500; % 覆蓋區域的長度,單位為公里
        ASDdeg = 10; % 角度分散,單位為度
        LEOPathloss_model = 'PathLoss'; % 設置為衛星路徑損耗模型
        DeltadB = 10; % 設置所需的 SNR 值(dB)

        %Compute channel statistics for one setup
        [R,channelGaindB] = functionExampleSetup(L,K,Mmax,accuracy,ASDdeg,squareLength,Pathloss_model,deployment,pathloss);
        
        %Compute the normalized average channel gain, where the normalization
        %is based on the noise power
        channelGainOverNoiseOriginal = channelGaindB - noiseVariancedBm;
        
        %Extract the average channel gains before power control
        channelGainOverNoise = channelGainOverNoiseOriginal;
        
        %Go through all cells
        for j = 1:L
            
            %Scale the average channel gains by applying the inverse statistical power control
            backoff = channelGainOverNoiseOriginal(:,j,j) + 10*log10(p) - DeltadB;
            channelGainOverNoise(:,j,:) = channelGainOverNoiseOriginal(:,j,:) - repmat(backoff,[1 1 L]);
            
        end
%這段程式碼用於遍歷所有細胞，並根據反向統計功率控制對平均通道增益進行縮放。
%對於每個細胞 j，首先從 channelGainOverNoiseOriginal 中提取平均通道增益的初始值。
%然後計算反向統計功率控制的偏移量 backoff。這是由原始通道增益 channelGainOverNoiseOriginal、功率 p 和 DeltadB（SNR 偏移）計算得到的。DeltadB 是 SNR 偏移（dB）的大小。
%將 backoff 應用於 channelGainOverNoiseOriginal，將其廣播到所有細胞上。這樣就可以獲得根據反向統計功率控制後的新的平均通道增益值，存儲在 channelGainOverNoise 中

        %Go through all number of antennas
        for m = 1:length(Mrange)
            
            %Output simulation progress
            disp([num2str(m) ' antennas out of ' num2str(length(Mrange))]);
            
            %Go through all pilot reuse factors
            for s = 1:length(fRange)
                
                %Extract pilot reuse factor
                f = fRange(s);
                
                if c == 1 %Correlated fading
                    
                    %Generate channel realizations with estimates and estimation
                    %error correlation matrices
                    disp('Running for correlated fading');
                    [Hhat,C,tau_p,Rscaled] = functionChannelEstimates(R(1:Mrange(m),1:Mrange(m),:,:,:),channelGainOverNoiseOriginal,DeltadB_pilot,nbrOfRealizations,Mrange(m),K,L,p,f,deployment);
                    %為什麼這邊要加DeltadB_pilot
                elseif c == 2 %Uncorrelated fading
                    
                    %Generate uncorrelated correlation matrices
                    Runcorr = repmat(eye(Mrange(m)),[1 1 K L L]);
%這行程式碼的意思是創建一個五維陣列 Runcorr，其內容是重複的 Mrange(m) x Mrange(m) 單位矩陣，
%並且在第三、第四和第五個維度上分別有 K、L 和 L 個這樣的單位矩陣 
                    %Generate channel realizations with estimates and estimation
                    %error correlation matrices
                    disp('Running for uncorrelated fading');
                    [Hhat,C,tau_p,Rscaled] = functionChannelEstimates(Runcorr,channelGainOverNoise,DeltadB_pilot,nbrOfRealizations,Mrange(m),K,L,p,f,deployment);
                    
                end

                %Compute SEs using Theorem 4.1
                [SE_MR,SE_RZF,SE_MMMSE,SE_ZF,SE_SMMSE] = functionComputeSE_UL(Hhat,C,Rscaled,tau_c,tau_p,nbrOfRealizations,Mrange(m),K,L,p);

                %Save average sum SE per cell
                sumSE_MR(ll,s,m,n) = mean(sum(SE_MR,1));
                sumSE_ZF(ll,s,m,n) = mean(sum(SE_ZF,1));
                sumSE_SMMSE(ll,s,m,n) = mean(sum(SE_SMMSE,1));
                sumSE_RZF(ll,s,m,n) = mean(sum(SE_RZF,1));
                sumSE_MMMSE(ll,s,m,n) = mean(sum(SE_MMMSE,1));
                
                %Delete large matrices
                clear Hhat C Rscaled;
                
            end
            
        end
        
        %Delete large matrices
        clear R;
        
    end
    
        save('../SimulationResults/Allen_simulation','deployment','pathloss','c','DeltadB','DeltadB_pilot','Pathloss_model','lambdaRange','fRange','Mrange','K','sumSE_MR','sumSE_ZF','sumSE_RZF','sumSE_SMMSE','sumSE_MMMSE','SE_MR','SE_RZF','SE_MMMSE','SE_ZF','SE_SMMSE');

end

 %% Plot the simulation results
 mm = length(Mrange);

 % Plot aggregate SE per cell vs BS density
for ff = 1:length(fRange)
    figure(length(lambdaRange)+ff);
    hold on; box on;
    plot(lambdaRange,squeeze(mean(sumSE_MMMSE(:,ff,mm,:),4)),'rd-','LineWidth',1,'DisplayName','MMMSE');
    plot(lambdaRange,squeeze(mean(sumSE_SMMSE(:,ff,mm,:),4)),'b:','LineWidth',1,'DisplayName','SMMSE');
    %plot(lambdaRange,squeeze(mean(sumSE_RZF(:,ff,mm,:),4)),'k-.','LineWidth',1,'DisplayName','RZF');
    plot(lambdaRange,squeeze(mean(sumSE_ZF(:,ff,mm,:),4)),'mo-','LineWidth',1,'DisplayName','ZF');
    plot(lambdaRange,squeeze(mean(sumSE_MR(:,ff,mm,:),4)),'bs-','LineWidth',1,'DisplayName','MR');
    legend('Location','best');
    xlabel('BS');
    ylabel('SE');
    title(['fRange = ' num2str(fRange(ff))]); % 設置標題，顯示當前 fRange 的值
end
 
