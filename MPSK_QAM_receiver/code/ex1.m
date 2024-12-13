%16QAM
%para
clear
fc = 3*10^(6);
Ts = 10^(-6);
Es = 10^(-6);
a = sqrt(Es); %wave amplitude
sample_freq = 10^(8);
sample_time = 1/sample_freq;
No_symbols = 1000;
t = [0:sample_time:Ts-sample_time];
phi1 = sqrt(2/Ts)*cos(2*pi*fc*t);
phi2 = sqrt(2/Ts)*sin(2*pi*fc*t); 
%generate m
select_1_3 = [1 3 -1 -3];
m = reshape(randsample(select_1_3, No_symbols*2, true), No_symbols, 2);
%generate R vactor (basis:phi1, phi2)
R = [];
SNR_db = 20;
SNR = 10.^(SNR_db/10); %Es/N0 = 20 (db)
N0 = Es/SNR;
for i = 1: No_symbols
    si_wave =  a*m(i,1)*phi1 + a*m(i,2) * phi2;
    siXphi1 = si_wave .* phi1; siXphi2 = si_wave .* phi2;
    int_siXphi1 = sample_time*sum(siXphi1);
    int_siXphi2 = sample_time*sum(siXphi2);
    %add noise (after correlator)
    Ri = [int_siXphi1 + normrnd(0,sqrt(N0/2)), int_siXphi2 + normrnd(0,sqrt(N0/2))];
    R = [R;Ri];
end
figure
scatter(R(:, 1), R(:, 2), 'x')
title(sprintf('16 QAM  SNR:%.2f dB', SNR_db))
xlabel('\phi_{1}(t)') 
ylabel('\phi_{2}(t)') 
grid on

%8 PSK
%para
clear
fc = 3*10^(6);
Ts = 10^(-6);
Es = 10^(-6);
a = sqrt(Es); %wave amplitude
M = 8; %QAM
sample_freq = 10^(8);
sample_time = 1/sample_freq;
No_symbols = 1000;
t = [0:sample_time:Ts-sample_time];
phi1 = sqrt(2/Ts)*cos(2*pi*fc*t);
phi2 = sqrt(2/Ts)*sin(2*pi*fc*t); 
%generate m
m = randsample([1:M], No_symbols, true);
%generate R vactor (basis:phi1, phi2)
R = [];
SNR_db = 20;
SNR = 10.^(SNR_db/10); %Es/N0 = 20 (db)
N0 = Es/SNR;
for i = 1: No_symbols
    si_wave =  sqrt(2*Es/Ts)*cos(2*pi*fc*t+ (m(i)-1)*2*pi/M);
    siXphi1 = si_wave .* phi1; siXphi2 = si_wave .* phi2;
    int_siXphi1 = sample_time*sum(siXphi1);
    int_siXphi2 = sample_time*sum(siXphi2);
    %add noise (after correlator)
    Ri = [int_siXphi1 + normrnd(0,sqrt(N0/2)), int_siXphi2 + normrnd(0,sqrt(N0/2))];
    R = [R;Ri];
end
figure
scatter(R(:, 1), R(:, 2), 'x')
title(sprintf('8 PSK  SNR:%.2f dB', SNR_db))
xlabel('\phi_{1}(t)') 
ylabel('\phi_{2}(t)') 
grid on






