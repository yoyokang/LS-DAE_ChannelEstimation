
function [Xtraining_Cell, Xtraining_Array, Ytraining_regression_array, Xvalidation_regression, Yvalidation_regression] =Link5_data_Generation(Training_set_ratio, SNRrange, Num_of_frame)

M = 4; % QPSK

Num_of_subcarriers = 63; %子载波数量
Num_of_FFT = Num_of_subcarriers + 1;%FFT点数
length_of_CP = 16;%循环前缀长度

Num_of_symbols = 1;%每帧包含的OFDM符号的个数
Num_of_pilot = 1;%每帧包含导频个数
Frame_size = Num_of_symbols + Num_of_pilot;%帧的大小,一"帧"包含一个OFDM符号和一个导频

Pilot_interval = Frame_size / Num_of_pilot;%导频间隔
Pilot_starting_location = 1;%导频起始位置

length_of_symbol = Num_of_FFT + length_of_CP;
Num_of_QPSK_symbols = Num_of_subcarriers * Num_of_symbols * Num_of_frame;
Num_of_bits = Num_of_QPSK_symbols * log2(M);

Xtraining_Cell = cell(Training_set_ratio * Num_of_frame, 1);
% Xtraining_Array = zeros(Num_of_FFT * Frame_size * 2, 1, 1, Training_set_ratio * Num_of_frame);
Xtraining_Array = zeros(Num_of_FFT * 2, 1, 1, Training_set_ratio * Num_of_frame);%num_of_FFT对应DAE输入神经元个数，实部虚部合起来为两倍
Ytraining_regression_array = zeros(1, 1, 128, Training_set_ratio * Num_of_frame);

% Xvalidation_regression = zeros(Num_of_FFT * Frame_size * 2, 1, 1, Num_of_frame - Training_set_ratio * Num_of_frame);
Xvalidation_regression = zeros(Num_of_FFT * 2, 1, 1, Num_of_frame - Training_set_ratio * Num_of_frame);
Yvalidation_regression = zeros(1, 1, 128, Num_of_frame - Training_set_ratio * Num_of_frame);

%winner2信道
frmLen = 160;

commSupportPackageCheck("CST_WINNER2");

AA(1) = winner2.AntennaArray("UCA",64,0.3);
AA(2) = winner2.AntennaArray("ULA",2,0.05);

BSIdx = {1}; % Index in antenna array inventory vector
MSIdx = [2];     % Index in antenna array inventory vector
numLinks = 1;            % Number of links
range = 10000;             % Layout range (meters)
seed = 101;
cfgLayout = winner2.layoutparset(MSIdx,BSIdx,numLinks,AA,range,seed);

cfgLayout.Pairing = [1 ; 2];  % Index in cfgLayout.Stations
cfgLayout.ScenarioVector = [10];     % 6 for B4, 11 for C2 and 13 for C4
cfgLayout.PropagConditionVector = [0];  % 0 for NLOS

numBSSect = sum(cfgLayout.NofSect);
numMS = length(MSIdx);

cfgLayout.Stations(1).Pos(1:2) = [5000; 1000];
cfgLayout.Stations(2).Pos(1:2) = [7000; 4000];

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

for Frame = 1:Num_of_frame

% 产生Bit流
N = Num_of_subcarriers * Num_of_symbols;
data = randi([0 1], N, 2);%一帧包含的bit数据
dataSym = bi2de(data);%二进制转为十进制，00-0，01-1，11-3，10-2，省去串并转换

% QPSK调制
N1 = size(dataSym, 1) * size(dataSym, 2);%计算dataSym中符号总数
QPSK_symbol = zeros(N1, 1);
QPSK_symbol(( dataSym == 0)) = 1 + 1j;%QPSK映射
QPSK_symbol(( dataSym == 1)) = -1 + 1j;
QPSK_symbol(( dataSym == 3)) = -1 - 1j;
QPSK_symbol(( dataSym == 2)) = 1 - 1j;

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
% chanInfo = info(winChannel);
% numTx = chanInfo.NumBSElements(1);
% Rs = chanInfo.SampleRate(1);

input_signal = repmat(Transmitted_signal, 1, 64);  % 复制 A，使其变成 1024x8 的矩阵，每列元素相同

y = winChannel(input_signal);

SNR = randi(SNRrange);
n = awgn(Transmitted_signal,SNR,'measured');
n = repmat(n, 1, 2);
x = n + y{1};

%分集接收
% 计算每根天线信号的幅度
mag_X11 = abs(x(:, 1));  % 计算第一列的幅度
mag_X12 = abs(x(:, 2));  % 计算第二列的幅度

% 计算最大比合并的加权系数
w1 = mag_X11.^2;  % 第一根天线的权重（幅度的平方）
w2 = mag_X12.^2;  % 第二根天线的权重（幅度的平方）

% 合并信号：通过加权求和
X1_total = (w1 .* y{1}(:, 1) + w2 .* y{1}(:, 2)) ./ (w1 + w2);  % 合并后的信号

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
% pilot_QPSK = pilot_QPSK.*(1-1j);
% pilot_datasym = ones(size(H_LS,1),size(H_LS,2)); 
% pilot_datasym = pilot_datasym.*2;

[feature_of_H_LS, label_of_regression] = Extract_Feature_DAE(H_LS, H_QPSK);

if Frame <= fix(Training_set_ratio * Num_of_frame)
    Xtraining_Cell{Frame, 1} = feature_of_H_LS;
    Xtraining_Array(:, 1, 1, Frame) = feature_of_H_LS;
    Ytraining_regression_array(1, 1, :, Frame) = label_of_regression;
else
    Xvalidation_regression(:, 1, 1, Frame - Training_set_ratio * Num_of_frame, 1) = feature_of_H_LS;
    Yvalidation_regression(1, 1, :, Frame - Training_set_ratio * Num_of_frame) = label_of_regression;
end

end

end
