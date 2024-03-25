function [PL,path_loss_dB]=PathLoss(fc,h_BS,d_2D)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Path loss model in urban micro NLOS scenario for mmwave
%--------------------------intput parameter--------------------------------
%   fc        : carrier frequency(GHz)
%   h_BS      : BS height[meter]
%   d_2D      : 2D distance between BS and UE(meter)
%--------------------------output parameter--------------------------------
%   PL        : Path loss(linear)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    h_UE=1.5+rand*21;                                                      % the height of UE(1.5m~21m)
    d_3D=norm([(h_BS-h_UE) d_2D],2);
    path_loss_dB=20*log10(d_3D)+log10(fc)+log10((4*3.14)/(3*10^8));% dB
    PL=10^(path_loss_dB/10);                                               % dB->linear
end