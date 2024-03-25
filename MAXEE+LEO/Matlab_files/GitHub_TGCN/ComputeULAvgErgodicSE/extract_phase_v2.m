function [exp_j_theta]=extract_phase_v2(Z,N)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% theta = angle(z) returns the phase angle in the interval [-£k,£k] for each element of a complex array z
% extract the phase of each element of matrix or vector to design our
% analog beamforming
% Z                 :input matrix or vector
% exp_j_theta       :output unit modulus matrix or vector
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    exp_j_theta=exp(1i * angle(Z))/sqrt(N);
end