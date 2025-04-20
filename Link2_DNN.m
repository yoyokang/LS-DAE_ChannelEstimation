SNR_Range = -20:20;%信噪比变化范围
%初始化矩阵
BER_over_SNR_LS = zeros(length(SNR_Range), 1);%用于存放不同信噪比下LS信道估计的误比特率
BER_over_SNR_MMSE = zeros(length(SNR_Range), 1);%用于存放不同信噪比下MMSE信道估计的误比特率
BER_over_SNR_DNN = zeros(length(SNR_Range), 1);
%

%winner2信道
frmLen = 160;

commSupportPackageCheck("CST_WINNER2");

AA(1) = winner2.AntennaArray("UCA",64,0.1);
AA(2) = winner2.AntennaArray("ULA",2,0.05);

BSIdx = {1}; % Index in antenna array inventory vector
MSIdx = [2];     % Index in antenna array inventory vector
numLinks = 1;            % Number of links
range = 10000;             % Layout range (meters)
seed = 101;
cfgLayout = winner2.layoutparset(MSIdx,BSIdx,numLinks,AA,range,seed);

cfgLayout.Pairing = [1;2];  % Index in cfgLayout.Stations
cfgLayout.ScenarioVector = [1];     % 6 for B4, 11 for C2 and 13 for C4
cfgLayout.PropagConditionVector = [0];  % 0 for NLOS

numBSSect = sum(cfgLayout.NofSect);
numMS = length(MSIdx);

% Set up positions for BS\MS sectors. 
cfgLayout.Stations(1).Pos(1:2) = [9500; 7000];
cfgLayout.Stations(2).Pos(1:2) = [9450; 7000];

for i = numBSSect + (1:numMS)
    cfgLayout.Stations(i).Velocity = rand(3,1) - 0.5;
end

cfgWim = winner2.wimparset;
cfgWim.NumTimeSamples = frmLen;
cfgWim.IntraClusterDsUsed = "yes";
cfgWim.CenterFrequency = 5.25e9;
cfgWim.UniformTimeSampling = "no";
cfgWim.ShadowingModelUsed = "yes";
cfgWim.PathLossModelUsed = "yes";
cfgWim.RandomSeed = 31415926;       % For repeatability

winChannel = comm.WINNER2Channel(cfgWim,cfgLayout);

chanInfo = info(winChannel);
numTx = chanInfo.NumBSElements(1);
Rs = chanInfo.SampleRate(1);

    x = [ones(1,numTx); zeros(frmLen-1,numTx)];
    h = winChannel(x);
    h1 = h{1}(:,1);  

for SNR = SNR_Range

M = 4; % QPSK
k = log2(M);

Num_of_subcarriers = 63; %子载波数量
Num_of_FFT = Num_of_subcarriers + 1;%FFT点数
length_of_CP = 16;%循环前缀长度

Num_of_symbols = 1;%每帧包含的OFDM符号的个数
Num_of_pilot = 1;%每帧包含导频个数
Frame_size = Num_of_symbols + Num_of_pilot;%帧的大小,一"帧"包含一个OFDM符号和一个导频
Num_of_frame = 10;%总共传输"帧"数

Pilot_interval = Frame_size / Num_of_pilot;%导频间隔
Pilot_starting_location = 1;%导频起始位置

length_of_symbol = Num_of_FFT + length_of_CP;%一个OFDM符号的长度

Num_of_QPSK_symbols = Num_of_subcarriers * Num_of_symbols * Num_of_frame;% 总共QPSK符号个数
Num_of_bits = Num_of_QPSK_symbols * log2(M);%总Bit数，用于计算总的误比特率

% Num_of_QPSK_symbols_DNN = 8 * Num_of_symbols * Num_of_frame;%DNN的QPSK符号个数
% Num_of_bits_DNN = Num_of_QPSK_symbols_DNN * k;%DNN总Bit数，用于计算总的误比特率

%初始化矩阵
numErrs_bit_LS = zeros(Num_of_frame, 1);%用于记录LS信道估计下每"帧"的误比特数
BER_in_frame_LS = zeros(Num_of_frame, 1);%用于记录LS信道估计下每"帧"的BER

numErrs_bit_MMSE = zeros(Num_of_frame, 1);%用于记录MMSE信道估计下每"帧"的误比特数
BER_in_frame_MMSE = zeros(Num_of_frame, 1);%用于记录MMSE信道估计下每"帧"的BER

numErrs_bit_DNN = zeros(Num_of_frame, 1);%用于记录MMSE信道估计下每"帧"的误比特数
BER_in_frame_DNN = zeros(Num_of_frame, 1);%用于记录MMSE信道估计下每"帧"的BER

Multipath_h = zeros(Num_of_frame * Frame_size * length_of_symbol, 1);

for Frame = 1:Num_of_frame

% 产生Bit流
N = Num_of_subcarriers * Num_of_symbols;
data = randi([0 1], N, 2);%一帧包含的bit数据
dataSym = bi2de(data);%DNN需要
QPSK_modulated = bi2de(data);%二进制转为十进制，00-0，01-1，11-3，10-2，省去串并转换
Data = reshape(data, [], 1);

% QPSK调制
N1 = size(QPSK_modulated, 1) * size(QPSK_modulated, 2);%计算dataSym中符号总数
QPSK_symbol = zeros(N1, 1);
QPSK_symbol(( QPSK_modulated == 0)) = 1 + 1j;%QPSK映射
QPSK_symbol(( QPSK_modulated == 1)) = -1 + 1j;
QPSK_symbol(( QPSK_modulated == 3)) = -1 - 1j;
QPSK_symbol(( QPSK_modulated == 2)) = 1 - 1j;
QPSK_signal = reshape(QPSK_symbol, Num_of_subcarriers, Num_of_symbols);

% 插入导频
Pilot_value = 1 - 1j;
Pilot_location = Pilot_starting_location : Pilot_interval : Frame_size;

data_location = 1 : Frame_size;
data_location(Pilot_location(:)) = [];

%初始化IFFT前存放数据的矩阵
data_in_IFFT = zeros(Num_of_FFT - 1, Frame_size);
%填入数据和导频
data_in_IFFT(:, Pilot_location(:)) = Pilot_value;
data_in_IFFT(:, data_location(:)) = QPSK_signal;
data_in_IFFT = [zeros(1, Frame_size); data_in_IFFT];

% 形成OFDM符号
data_in_CP = ifft(data_in_IFFT);% IFFT
data_in_CP = sqrt(Num_of_FFT) * data_in_CP;
% 加CP
Signal_from_baseband = [data_in_CP(Num_of_FFT - length_of_CP + 1 : end, :); data_in_CP];
Transmitted_signal = reshape(Signal_from_baseband, [], 1);

input_signal = repmat(Transmitted_signal, 1, 64);  % 复制 A，使其变成 1024x8 的矩阵，每列元素相同
y = winChannel(input_signal);

% n = awgn(y{1},SNR,'measured');
n = awgn(Transmitted_signal,SNR,'measured');
n = repmat(n, 1, 2);
y{1}= n + y{1};
x = y{1};

%分集接收
% 计算每根天线信号的幅度
mag_X11 = abs(x(:, 1));  % 计算第一列的幅度
mag_X12 = abs(x(:, 2));  % 计算第二列的幅度

% 计算最大比合并的加权系数
w1 = mag_X11.^2;  % 第一根天线的权重（幅度的平方）
w2 = mag_X12.^2;  % 第二根天线的权重（幅度的平方）

% 合并信号：通过加权求和
X1_total = (w1 .* x(:, 1) + w2 .* x(:, 2)) ./ (w1 + w2);  % 合并后的信号

% 去除CP
Received_Signal = reshape(X1_total, [], Frame_size);
Received_Signal_CP_removed = Received_Signal(length_of_CP + 1 : end, :);

% FFT变换
Received_Signal_FFT = fft(Received_Signal_CP_removed) / sqrt(Num_of_FFT);

% 提取导频符号
Received_Pilot = Received_Signal_FFT(:, Pilot_location(:));

% 已知的导频值
Pilot_value = 1 - 1j;

% LS信道估计
H_LS = Received_Pilot ./ Pilot_value;
H_QPSK = ones(size(H_LS,1),size(H_LS,2)); 

% MMSE信道估计
% h1 = h{1}(:,1);

Nfft = Num_of_FFT;
H_MMSE_h = zeros(Nfft, Frame_size);
SNR_HEX = 10^(SNR / 10);
Np = Nfft;
% H_LS = Received_Pilot ./ Pilot_value;
Nps = 1;
h = h1';
k = 0: length(h) - 1;
hh = h * h';
tmp = h .* conj(h) .* k;
r = sum(tmp) / hh;
r2 = tmp * k .'/hh;

tau_rms = sqrt(r2 - r^2);
df = 1/Nfft;
j2pi_tau_df = 1j * 2 * pi * tau_rms * df;
K1 = repmat([0 : Nfft - 1].', 1, Np);
K2 = repmat([0 : Np - 1], Nfft, 1);
rf = 1./(1 + j2pi_tau_df * (K1 - K2 * Nps));
K3 = repmat([0 : Np - 1].', 1, Np);
K4 = repmat([0 : Np - 1], Np, 1);
rf2 = 1./(1 + j2pi_tau_df * Nps * (K3 - K4));
Rhp = rf;
Rpp = rf2 + (eye(length(H_LS)) / SNR_HEX);
H_MMSE = Rhp * pinv(Rpp) * H_LS;

for i_MMSE = 1 : Frame_size
    H_MMSE_h(:, i_MMSE) = H_MMSE;
end

H_MMSE = H_MMSE_h;

%DNN处理数据
%  load('DNN_Trained.mat');
load('DAE_Link2.mat');
% 只取8个QPSK符号
[feature_of_H_LS, ~] = Extract_Feature_DAE(H_LS, H_QPSK);
Received_H = predict(DNN_Trained, feature_of_H_LS);
Received_H = transpose(Received_H);
H_DAE = Received_H(1:2:end, :) + 1j * Received_H(2:2:end, :);

%传统方法信道均衡
Received_Signal_LS = Received_Signal_FFT./ H_LS;
Received_Signal_MMSE = Received_Signal_FFT ./ H_MMSE;
%加入DAE
Received_Signal_DAE = Received_Signal_FFT./ H_DAE;

% 提取数据符号
Received_Signal_L= Received_Signal_LS(2:end, data_location(:));
Received_Signal_M= Received_Signal_MMSE(2:end, data_location(:));
Received_Signal_D= Received_Signal_DAE(2:end, data_location(:));

% QPSK解调
N2 = size(Received_Signal_L, 1) * size(Received_Signal_L, 2);
N3 = size(Received_Signal_M, 1) * size(Received_Signal_M, 2);
N4 = size(Received_Signal_D, 1) * size(Received_Signal_D, 2);

QPSK_demodulated_LS = zeros(N2, 1);
QPSK_demodulated_LS(find((angle(Received_Signal_L) > 0) .* (angle(Received_Signal_L) < (pi * 1/2)))) = 0;
QPSK_demodulated_LS(find((angle(Received_Signal_L) > (pi/2)) .* (angle(Received_Signal_L) < pi))) = 1;
QPSK_demodulated_LS(find((angle(Received_Signal_L) > (- pi)) .* (angle(Received_Signal_L) < (- pi/2)))) = 3;
QPSK_demodulated_LS(find((angle(Received_Signal_L) > (- pi/2)) .* (angle(Received_Signal_L) < 0))) = 2;

QPSK_demodulated_MMSE = zeros(N3, 1);
QPSK_demodulated_MMSE(find((angle(Received_Signal_M) > 0) .* (angle(Received_Signal_M) < (pi * 1/2)))) = 0;
QPSK_demodulated_MMSE(find((angle(Received_Signal_M) > (pi/2)) .* (angle(Received_Signal_M) < pi))) = 1;
QPSK_demodulated_MMSE(find((angle(Received_Signal_M) > (- pi)) .* (angle(Received_Signal_M) < (- pi/2)))) = 3;
QPSK_demodulated_MMSE(find((angle(Received_Signal_M) > (- pi/2)) .* (angle(Received_Signal_M) < 0))) = 2;

QPSK_demodulated_DNN = zeros(N4, 1);
QPSK_demodulated_DNN(find((angle(Received_Signal_D) > 0) .* (angle(Received_Signal_D) < (pi * 1/2)))) = 0;
QPSK_demodulated_DNN(find((angle(Received_Signal_D) > (pi/2)) .* (angle(Received_Signal_D) < pi))) = 1;
QPSK_demodulated_DNN(find((angle(Received_Signal_D) > (- pi)) .* (angle(Received_Signal_D) < (- pi/2)))) = 3;
QPSK_demodulated_DNN(find((angle(Received_Signal_D) > (- pi/2)) .* (angle(Received_Signal_D) < 0))) = 2;

% 将QPSK符号转换为比特流
data_demodulated_LS = de2bi(QPSK_demodulated_LS, 2);
Data_demodulated_LS = reshape(data_demodulated_LS, [], 1);
data_demodulated_MMSE = de2bi(QPSK_demodulated_MMSE, 2);
Data_demodulated_MMSE = reshape(data_demodulated_MMSE, [], 1);
data_demodulated_DNN = de2bi(QPSK_demodulated_DNN, 2);
Data_demodulated_DNN = reshape(data_demodulated_DNN, [], 1);
%计算一帧的BER
numErrs_bit_LS(Frame, 1) = sum(sum(round(Data) ~= round(Data_demodulated_LS)));%存储一帧的误比特数
BER_in_frame_LS(Frame, 1) = numErrs_bit_LS(Frame, 1) / length(Data);%存储一帧的BER

numErrs_bit_MMSE(Frame, 1) = sum(sum(round(Data) ~= round(Data_demodulated_MMSE)));%存储一帧的误比特数
BER_in_frame_MMSE(Frame, 1) = numErrs_bit_MMSE(Frame, 1) / length(Data);%存储一帧的BER

numErrs_bit_DNN(Frame, 1) = sum(sum(round(Data) ~= round(Data_demodulated_DNN)));%存储一帧的误比特数
BER_in_frame_DNN(Frame, 1) = numErrs_bit_DNN(Frame, 1) / length(Data);%存储一帧的BER

end

%计算所有帧的BER
BER_LS = sum(numErrs_bit_LS, 1) / Num_of_bits;
BER_MMSE = sum(numErrs_bit_MMSE, 1) / Num_of_bits;
BER_DNN = sum(numErrs_bit_DNN, 1) / Num_of_bits;

BER_over_SNR_LS(SNR+21, 1) = BER_LS;%存储BER数据
BER_over_SNR_MMSE(SNR+21, 1) = BER_MMSE;
BER_over_SNR_DNN(SNR+21, 1) = BER_DNN;

end

plot(SNR_Range, BER_over_SNR_LS, 'b-o', SNR_Range, BER_over_SNR_MMSE, 'r-*',SNR_Range, BER_over_SNR_DNN, 'g-^');
xlabel('SNR (dB)');
ylabel('BER');
legend('LS', 'MMSE','LS+DAE');
title("WINNER信道A1")
grid on;


