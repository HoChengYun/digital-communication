%16QAM
%para
clear
fc = 3*10^(6);
Ts = 10^(-6);
Es = 10^(-6);
a = sqrt(Es); %wave amplitude
M = 16;
sample_freq = 10^(8);
sample_time = 1/sample_freq;
No_symbols = 100000;
t = [0:sample_time:Ts-sample_time];
phi1 = sqrt(2/Ts)*cos(2*pi*fc*t);
phi2 = sqrt(2/Ts)*sin(2*pi*fc*t); 
%generate m
m = randsample([1:M], No_symbols, true);
QAM_conj = [1 1;1 -1; -1 1;-1 -1;...
         3 1;3 -1; -3 1;-3 -1;...
         1 3;1 -3; -1 3;-1 -3;...
         3 3;3 -3; -3 3;-3 -3];
m_conj = QAM_conj(m, :);
%set snr
SNR_db = [-6:2:20];
SNR = 10.^(SNR_db/10); %Es/N0 = 20 (db)
N0_array = Es./SNR;
%bit error rate
Pe = [];
%generate S vector (basis:phi1, phi2)
S = [];
for i = 1: M
    si_wave =  a*QAM_conj(i,1)*phi1 + a*QAM_conj(i,2)*phi2;
    siXphi1 = si_wave .* phi1; siXphi2 = si_wave .* phi2;
    int_siXphi1 = sample_time*sum(siXphi1);
    int_siXphi2 = sample_time*sum(siXphi2);
    Si = [int_siXphi1, int_siXphi2];
    S = [S ;Si];
end
%generate R vactor (basis:phi1, phi2)
for N0 = N0_array
    R = [];
    D = zeros(No_symbols,M);
    D_min = [];
    m_hat = [];
    for i = 1: No_symbols
        si_wave =  a*m_conj(i,1)*phi1 + a*m_conj(i,2)*phi2;
        siXphi1 = si_wave .* phi1; siXphi2 = si_wave .* phi2;
        int_siXphi1 = sample_time*sum(siXphi1);
        int_siXphi2 = sample_time*sum(siXphi2);
        %add noise (after correlator)
        Ri = [int_siXphi1 + normrnd(0, sqrt(N0/2)), int_siXphi2 + normrnd(0,sqrt(N0/2))];
        R = [R;Ri];
    end
    %generate Di
    for i = 1:M
        distance = -abs((R-S(i,:)).^(2));
        D(:, i) = distance(:, 1) + distance(:, 2);
    end
    D = D.';
    %find m_hat
    [D_max, m_hat] = max(D);
    %ber generate
    bit_errorNUM = sum((m - m_hat) ~= 0);
    BER = bit_errorNUM/No_symbols;
    Pe = [Pe BER];
end
%calculate error rate theory
Pe_theory = 3*qfunc(sqrt(4*a*a./(2*N0_array))) - (9/4)*qfunc(sqrt(4*a*a./(2*N0_array))).^(2);
figure
semilogy(SNR_db, Pe, '-x', SNR_db, Pe_theory, '-o')
title('bit error rate 16QAM'); 
xlabel('Es/N0 (db)'); 
ylabel('BER'); 
legend('experiment','theory','Location', 'southwest');