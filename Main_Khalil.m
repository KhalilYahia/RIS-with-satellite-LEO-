%% This is the main file
% This created by KHALIL YAHIA 
% khalilmohyahia@gmail.com
% MTUCI university, Moscow - Rusia
% Tishreen university, Latakia - Syria
% 29/02/2024
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
clear all;
fc = 2e9;
c = physconst('lightspeed');
lambda = c/fc;
rng(2023);
NoiseSD=-203; % Noise power spectral density db/Hz
BW = 20 *1e6; % Band width 20 MHz
N0 = NoiseSD+10*log10(BW); % dbW
Pt=10; % Transmitted power dbW
transmitted_Power = 10^(Pt/10); % transmitted Power in W
%t=[-523.97760176 -407.5381347  -291.09866765 -174.65920059  -58.21973353    0.           58.21973353  174.65920059  291.09866765  407.5381347 523.97760176];
[max_d t0]=Max_d_as_elevation_khalil(20,1500*1e3);
t=-t0:20:t0;
% Setup surface
M = 50;
N = 50;
% M_dx=2;
% N_dy=2;
dx =0.2*lambda;
dy = 0.2*lambda;

% Rotation around the x axis at an angle of 45 degrees
root_=45;
%% Create RIS with coordinates (0,0,0) 
RIS = RIS_Object(M,N,dx,dy,0,0,0,root_);

%% Create satellite
% Calculate satellite position using this function [Px Py Pz]= Calculate_SatePos_Khalil(t,rot)
[Sate_Location_x Sate_Location_y Sate_Location_z] = Calculate_SatePos_Khalil(0,0);
Sate_pos= [Sate_Location_x;Sate_Location_y;Sate_Location_z];

%% Create user
% User position (10,0,500) m
[Rec_pos_x Rec_pos_y Rec_pos_z] = Calculate_UserPos_Khalil(0,0);
User_pos = [Rec_pos_x;Rec_pos_y;Rec_pos_z];

%% A direct

SNR_withPlane=[];
SNR_withOutPlane=[];
theta_=[];
Gain_Tr=[];
channel_cap=[];
msg_lengh=1e4;
% 
ElevationAngle_user = Get_Elevation(Sate_pos,User_pos,1500*1e3);
N0 = N0*ones(1,msg_lengh);
for index_t=1:length(t)
    % Generate Rician channel
    h = Get_RicianF_Cha(fc,BW,ElevationAngle_user,msg_lengh);
    % Generate the direct channel
    [a tao R_direct]=ampgain_tao_direct(lambda,t(index_t),root_,transmitted_Power);
    R_direct = h .* (R_direct);
    noise_=(10.^(N0/10));
    PowerOfH_direct= abs(R_direct).^2;
    SNR_withOutPlane(index_t,:) = PowerOfH_direct ./noise_;
    % The Channel through RIS
    [A_ris tao_ris R_ris_and_Direct]=RIS.ampgain_tao_Ris_allelements(lambda,t(index_t),root_,tao,transmitted_Power,R_direct);
   
   
    PowerOfH_RIS_And_Direct = abs(R_ris_and_Direct + R_direct ).^2;
    SNR_withPlane(index_t,:) = (PowerOfH_RIS_And_Direct) ./noise_;
    
   
%     channel_cap(index_t)=10*log2(1+SNR_withPlane(index_t));

    %%%%%%%%
%     [Sate_Location_x Sate_Location_y Sate_Location_z] = Calculate_SatePos_Khalil(t(index_t),root_);
%     Sate_pos= [Sate_Location_x;Sate_Location_y;Sate_Location_z];
%     theta_(index_t)=theta_calc(Sate_pos,User_pos);
%     Gain_Tr(index_t)=Trans_Gain(theta_(index_t));
    
end

SNR_withOutPlane= 10*log10(mean(SNR_withOutPlane,2)); 
SNR_withPlane= 10*log10(mean(SNR_withPlane,2)); 

% 
%% 
p1=plot(t,SNR_withPlane,"Color","r",'LineWidth', 1.5,'DisplayName','SNR с использованием (RIS) 50x50 элементов');
hold on;
% p2=plot(t,SNR_withOutPlane,"Color","b",'LineWidth', 1.5,'DisplayName','SNR Без (RIS)');
ylabel('SNR (ДБ)')
xlabel('Время (секунда)')
legend

% figure
% 
% plot(t,channel_cap);
