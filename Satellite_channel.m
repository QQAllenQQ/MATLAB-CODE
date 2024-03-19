function [H_user] =Satellite_channel(Nt_user,Nr_BS,no_of_subcarriers,no_of_users,Nray,angspread,path_loss_dB,fc)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%--------------------------intput parameter--------------------------------
 Nt_user = 100            %Number of tx antenna per user
 Nr_BS = 128              %Number of rx antenna at BS
 no_of_subcarriers = 16   %Number of subcarriers
 no_of_users = 6         %Number of users
 Nray = 10               %Number of scatter(Arrived at the same time)
 angspread = pi/180*10   %angular spread(radian)
%--------------------------output parameter--------------------------------
%H_user             : %Frequency domain channel at kth subcarrier for all user
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Nscatter = Nray*Nray;                                                       % Total scatter
%% Initialzation matrix
txang = zeros(1,Nscatter);                                                 % Ray angle of tx
rxang = zeros(1,Nscatter);                                                 % Ray angle of rx
txarray = zeros(Nt_user,Nscatter);                                         % Array response vector of tx
rxarray = zeros(Nr_BS,Nscatter);                                           % Array response vector of rx
h_ray = complex(zeros(Nr_BS,Nt_user,Nray));                                % Time domain channel for each scatter for one user
h_ay = complex(zeros(Nr_BS,Nt_user,Nray));                                  
h_ay_DFT = complex(zeros(Nr_BS,Nt_user,Nray));                              % kth subcarrier frequency projection for one user (By DFT)
Hk = complex(zeros(Nr_BS,Nt_user,no_of_subcarriers));                      % Frequency domain channel at kth subcarrier for one user
H_user = complex(zeros(Nr_BS,Nt_user,no_of_subcarriers,no_of_users));      % Frequency domain channel at kth subcarrier for all user
if Nt_user <= 0 || Nr_BS <= 0 || no_of_subcarriers <= 0 || no_of_users <= 0 || Nray <= 0
    error('Input parameters must be positive integers.')
end
%% Channel model
for nu=1:no_of_users
    %--------------------------------------------------------------------------------------------------------------------------%   
    tx_ang_azimuth = unifrnd(0,2*pi,1,Nt_user);                                % Azimuth angle of transmitter (phi_t_mean )     , (uniform distribution) --> ULA
    tx_ang_elevation = unifrnd(0,pi,1,Nt_user);                              % Elevation angle of transmitter (theta_t_mean ) , (uniform distribution) --> UPA
    rx_ang_azimuth = unifrnd(0,2*pi,1,Nr_BS);                                % Azimuth angle of receiver (phi_r_mean )        , (uniform distribution) --> ULA
    rx_ang_elevation = unifrnd(0,pi,1,Nr_BS);                              % Elevation angle of receiver (theta_r_mean )    , (uniform distribution) --> UPA
    disp(['Processing user ', num2str(nu), '...']);
    %------------------------------------------------------------- Scatter angle -------------------------------------------------------------%   
   for t = 1:Nray
    if (t-1)*Nray+1 <= Nscatter
        tsc = max((t-1)*Nray+1, 1):min(t*Nray, Nscatter);
        tsc(tsc > Nscatter | tsc < 1) = [];
    else
        tsc = 1:Nscatter;
    end
    if  max(tsc) <= size(tx_ang_azimuth, 2) && max(tsc) <= size(tx_ang_elevation, 2)
        txang(1,tsc) = laprnd(1,Nray,tx_ang_azimuth(1,t),angspread);     % Aagle of departure (phi_t )   , (Laplacian distribution) --> ULA
        txang(2,tsc) = laprnd(1,Nray,tx_ang_elevation(1,t),angspread);    % Aagle of departure (theta_t ) , (Laplacian distribution) --> UPA
        rxang(1,tsc) = laprnd(1,Nray,rx_ang_azimuth(1,t),angspread);      % Aagle of arrival (phi_r)     , (Laplacian distribution) --> ULA
        rxang(2,tsc) = laprnd(1,Nray,rx_ang_elevation(1,t),angspread);    % Aagle of arrival (theta_r )   , (Laplacian distribution) --> UPA
    else
        disp(['For t = ', num2str(t), ', tsc = ', num2str(tsc)]);
        error('Index exceeds matrix dimensions.')
    end
end
    %-------------------------------------------------------- Array response vector ---------------------------------------------------------%
    ta = 1;                                                                                                                    % UPA ... 
       ra = 1;                                                                                                                    % 
       for w = 1:Nscatter                                                                                                         %
           for xx = 0:(sqrt(Nt_user)-1)                                                                                                 %
               for yy = 0:(sqrt(Nt_user)-1)                                                                                             %  
                   txarray(ta,w) = exp(-1i*pi*(xx*sin(txang(1,w))*cos(txang(2,w))+yy*cos(txang(2,w))))/sqrt(Nt_user);                     %        
                   ta = ta+1;                                                                                                     %
               end                                                                                                                %
           end                                                                                                                    %
        ta = 1;                                                                                                                   %
       end                                                                                                                        %
        for w = 1:Nscatter                                                                                                        %
            for xxx = 0:(sqrt(Nr_BS)-1)                                                                                           %   
                for yyy = 0:(sqrt(Nr_BS)-1)                                                                                       % 
                   rxarray(ra,w) = exp(-1i*pi*(xxx*sin(rxang(1,w))*cos(rxang(2,w))+yyy*cos(rxang(2,w))))/sqrt(Nr_BS);                %
                ra = ra+1;                                                                                                        %
                end                                                                                                               %
            end                                                                                                                   %
        ra = 1;                                                                                                                   %
        end
 end
    %% ------------------------------------------------------------------- channel model   ------------------------------------------------------------------- %    
    DS = 66e-9;
    r_tau = 2.1;                                            % Delay scaling parameter 
    Xn = rand(1,Nray);                                       % U~[0,1)
    tau_p = -r_tau*DS*log(Xn);                              % Exponential delay distribution
    tau = sort(tau_p-min(tau_p));                           % Normalize the delays by subtracting the minimum delay and sort them in ascending order 
    R = 500e3;            % Orbital altitude of LEO satellite
    h = 10;                % User height
    r = R + h;            % User-satellite distance
    d_LOS = r;
    fc = 28e9; % 28 GHz
    c = physconst('LightSpeed');     % speed of light (m/s)
    lambda = c/fc;      % wavelength (m)
    Kr = 10;
    GTX = 7; %dB
    GRX = 0; %dB
    LPat = 0.017; %dB
    path_loss_dB = 10*log10(10);
    Alpha = GTX+GRX-LPat-path_loss_dB;
    gamma = sqrt(1/Nray);
    u = sqrt(Alpha*Kr/(2*(Kr+1)));%mean
    d = sqrt(Alpha/(2*(Kr+1)));%variance
    ray_gain = gamma*[normrnd(u,d,1,Nscatter) + 1i*normrnd(u,d,1,Nscatter)]/sqrt(2);                          % Complex gain of each scatter
    ray_gain2 = [normrnd(0,1,1,Nscatter) + 1i*normrnd(0,1,1,Nscatter)]/sqrt(2)+(normrnd(sqrt(Kr),1,1,Nscatter));
    disp('Sample values of ray_gain:');
    disp(ray_gain(1:5));
    disp('Sample values of ray_gain2:');
    disp(ray_gain2(1:5));    
    %-------------------------------------------------------------------------%  
       max_i = min([size(ray_gain2, 1), size(rxarray, 1), size(txarray, 1)]);

for ray = 1:Nray
           scat = ray;                                                                         % Scatter Index
           temp = (sqrt(Alpha*Kr/(Kr+1)))*ray_gain2(:, scat)*exp(-1j*2*pi*d_LOS/lambda)*rxarray(:, 1)*txarray(:, 1)'*exp(-1j*2*pi*fc*tau(1)) + (sqrt(Alpha/(Kr+1)))*ray_gain(scat)*rxarray(:, scat)*txarray(:, scat).'*exp(-1j*2*pi*fc*tau(ray));
           h_ray(:, :, ray) = temp;
           disp('Size of h_ray:');
           disp(size(h_ray));
           disp('Sample values from h_ray:');
           disp(h_ray(:, :, 1));
end
    h_ay = sum(h_ray, 4);                                                                                     
        
    for k = 1: no_of_subcarriers
       for ray = 1: Nray
          h_ay_DFT(:,:,ray) = h_ay(:,:,ray)*exp(-1i*2*pi*(ray-1)*(k-1)/no_of_subcarriers);                    % kth subcarrier frequency projection for one user (By DFT)
       end                           
       Hk(:,:,k) = reshape(sum(h_ay_DFT,3), [Nr_BS, Nt_user]);                                                                           % Frequency domain channel at kth subcarrier for one user
    end
       for nu = 1:no_of_users
        H_user(:,:,:,nu) = Hk; % Assigning Hk to the nu-th user              % Frequency domain channel at kth subcarrier for all user
    end 
end % end of user index
