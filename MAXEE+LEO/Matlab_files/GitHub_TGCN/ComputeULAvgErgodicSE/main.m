%This Matlab script is the main of the article:
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

% �w�q satellite_channel ��ƪ��Ѽ�
Nt_user = 100;            % �C�ӥΤ᪺�o�g�ѽu�ƶq
Nr_BS = 128;              % �򯸪������ѽu�ƶq
no_of_subcarriers = 16;   % �l���i���ƶq
no_of_users = 6;          % �Τ᪺�ƶq
Nray = 10;                % �g�u�δ��g�骺�ƶq�]�P�ɨ�F�^
angspread = pi/180*10;    % ���׮i�}�]�H���׬����^
path_loss_dB = 100;       % ���|�l�ӡ]�ܨҤ����ȡ^
fc = 2e9;                 % ���i�W�v�]�Ҧp�A2 GHz�^

% �ϥΩw�q���Ѽƽե� satellite_channel ���
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

%�������W�]�G�������^�q�򯸡]BS�^��Τ�]�ơ]UE�^�������Z��(KM)
d_2D = 1000;

%���a���]BS�^������(Meter)
h_BS = 10;

%�b�o�q�{���X���Adeployment �ܼƥΨӨM�w��a�x�����p�覡�A����ؿ�ܡG
%�T�w����]Fixed Grid�^�G
%�� deployment �ܼƪ��Ȭ� 1 �ɡA��ܰ򩳯��b�@�өT�w������W�i�泡�p�C �o�N���۰�a�x����m�Q�w���T�w�A�åB�q�`�H�W�h�����檬�覡���G�b�ҿ�ϰ줺�C
%���êy�Q�I�L�{�]Homogeneous Poisson Point Process�AH-PPP�^�G

%�� deployment �ܼƪ��Ȭ� 2 �ɡA��ܰ�a�x�̷ӧ��êy�Q�I�{�Ƕi�泡�p�C �b�o�س��p�覡�U�A��a�x����m�O�H�����A�åB��`�y�Q�I�L�{���έp�S�ʡA�Y�b���w�ϰ줺�����N��m����a�x�ƶq�A�q�y�Q���G�C
%�ھڼ������ݭn�M��ڱ��p�A�i�H��ܦX�A����a�x���p�覡�C �T�w����A�Ω�ݭn�����a�x��m�è㦳�W�ߩʪ������A��H-PPP�A�Ω��u��M�H���������C


% Generates the large scale path-loss coefficient using (1) single slope or
% (2) multi slope model
pathloss = 1;

% (1) Correlated or (2) Uncorellated channel fading
c = 1;

% Density of BSs for the random deployment
lambdaRange = fliplr([1 2:2:10]);

%�o�q�{���X�O�Ψөw�q�H�����p���򯸡]BS�^���K�סC�b�o�ӱ��ҤU�A`lambdaRange` �O�@�ӦV�q�A���]�t�F���P���򯸱K�׭ȡA�o�ǭȷ|�b��u���ϥΡC
%�b�o�ӯS�w���Ҥl���A`fliplr([1 2:2:8])` ���N��O�إߤ@�ӦV�q�A�ӦV�q�]�t�F�q 1 �}�l���_�ƼƦC `[1 3 5 ...]`�A�M��ϥ� `fliplr` �禡�N�o�ӦV�q����A�ҥH�̫᪺ `lambdaRange` �N�ܦ��F `[1 1 1]`�C
%�o��ܦb�o�ӥ�u���A�򯸱K�׬O�T�w���B�� 1�A�Y�C���褽�����@�Ӱ򯸡C�o�˪��]�m�i��Ω�b�K���������U�i��į�����A�䤤�C�Ӱ򯸪��A�Ƚd��Q�X�z�a���|�A�H���ѥR�����л\�M�e�q�C

% Check the lengths of the variables
disp(length(lambdaRange));
disp(length(fRange));
disp(length(Mrange));

%Load fixed Propagation parameters %�o�q�O�I�sSetPropagationParameters�o���{���X
SetPropagationParameters;
PathLoss(fc,h_BS,d_2D);
%%

%Prepare to save simulation results
sumSE_MR = zeros(length(lambdaRange),length(fRange),length(Mrange),nbrOfSetups);
sumSE_ZF = zeros(length(lambdaRange),length(fRange),length(Mrange),nbrOfSetups);
sumSE_RZF = zeros(length(lambdaRange),length(fRange),length(Mrange),nbrOfSetups);
sumSE_SMMSE = zeros(length(lambdaRange),length(fRange),length(Mrange),nbrOfSetups);
sumSE_MMMSE = zeros(length(lambdaRange),length(fRange),length(Mrange),nbrOfSetups);

%�o�q�{���X�O�s��Allen_simulation��
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
    
% �o�q MATLAB �ϥΤF�@�� for �`���ӹM�� lambdaRange �V�q,�ÿ�X��e�������i�סC���\�ڬ��z�ԲӸ���:
%1. for ll = 1:length(lambdaRange) �w�q�F�@�Ӵ`��,�ܶq ll �q 1 ���N�� lambdaRange �����סC
%2. �b�C�����N��,ll ���� lambdaRange �����@�Ӱ򯸱K�׭ȡC
%3. disp() �Ω�b MATLAB �R�O������X�H���C
%4. num2str() ��ƥi�H��Ʀr�ഫ���r�Ŧ�C
%5. \[\] �Ω�r�Ŧ�s��,��{�����r�Ŧ�榡�ƪ��ĪG�C
%6. ��X���r�Ŧ�N���e�O�� ll �� lambda ��,�`�@�� length(lambdaRange) �ӡC
%����N��N�O,�b�`�����N lambdaRange ���C�@�Ӱ򯸱K�׮�,��X��e�������i�סC�o�˥i�H��K�l�ܪ��ɶ��B�檺�������窺�i�i���p�C
%�o�� Simple ���i�׿�X�b���ɶ��B�檺���礤�O�D�`���Ϊ�,�J�i�H�d�ݷ�e���A,�]�i�H���p�ѧE�ɶ��C�`��W�ܦn�a��{�F MATLAB �b�u�{���ΤW��X�H���M�ǹF���A���K�Q�ʡC
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
%�o�q�{���X�O�@�Ӱj��A�Ω�M�������������P�]�w�C�b�C�ӳ]�w���A���|�ˬd���p�������O�_�� 1 �� 2�C�ھڳ��p�������A���|�]�m�Ѽ� L ���ȡC
%�p�G���p�����O 1�A�h L �Q�]�m�� 16�C�o�N���ۦb�o�ر��p�U�A�t�Τ��� 16 �Ӱ򯸡]�βӭM�^�C
%�p�G���p�����O 2�A�h���ϥ� Poisson ���G�ͦ��@���H���� BS �ƶq L�A�o�Ӽƶq�O�ھڬY�ӽd�򤺪����� BS �ƱK�שM��ΰϰ쪺�j�p�p��o�X���C�ͦ��� L �ȥi�ण�O��ơA�]���ϥ� round ��ƶi��|�ˤ��J�C�P�ɡA���T�O L ���|���s�A�p�G�ͦ����Ȭ��s�A�h���s�ͦ��C
%�`���A�o�q�{���X�T�O�b�C�Ӽ����]�w���A�򯸡]�βӭM�^���ƶq L ���|�ھڳ��p�����ӽT�w�A�o�O���������@�ӭ��n�ѼơC


        squareLength = 500; % �л\�ϰ쪺����,��쬰����
        ASDdeg = 10; % ���פ���,��쬰��
        LEOPathloss_model = 'PathLoss'; % �]�m���ìP���|�l�Ӽҫ�
        DeltadB = 10; % �]�m�һݪ� SNR ��(dB)

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
%�o�q�{���X�Ω�M���Ҧ��ӭM�A�îھڤϦV�έp�\�v����省���q�D�W�q�i���Y��C
%���C�ӲӭM j�A�����q channelGainOverNoiseOriginal �����������q�D�W�q����l�ȡC
%�M��p��ϦV�έp�\�v��������q backoff�C�o�O�ѭ�l�q�D�W�q channelGainOverNoiseOriginal�B�\�v p �M DeltadB�]SNR �����^�p��o�쪺�CDeltadB �O SNR �����]dB�^���j�p�C
%�N backoff ���Ω� channelGainOverNoiseOriginal�A�N��s����Ҧ��ӭM�W�C�o�˴N�i�H��o�ھڤϦV�έp�\�v����᪺�s�������q�D�W�q�ȡA�s�x�b channelGainOverNoise ��

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
                    %������o��n�[DeltadB_pilot
                elseif c == 2 %Uncorrelated fading
                    
                    %Generate uncorrelated correlation matrices
                    Runcorr = repmat(eye(Mrange(m)),[1 1 K L L]);
%�o��{���X���N��O�Ыؤ@�Ӥ����}�C Runcorr�A�䤺�e�O���ƪ� Mrange(m) x Mrange(m) ���x�}�A
%�åB�b�ĤT�B�ĥ|�M�Ĥ��Ӻ��פW���O�� K�BL �M L �ӳo�˪����x�} 
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
    title(['fRange = ' num2str(fRange(ff))]); % �]�m���D�A��ܷ�e fRange ����
end
 
