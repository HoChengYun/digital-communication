%QPSK
%para
clear
fc = 3*10^(6);
Ts = 10^(-6);
Es = 1;
a = sqrt(Es); %wave amplitude
M = 4;
sample_freq = 10^(8);
sample_time = 1/sample_freq;
t = [0:sample_time:Ts-sample_time];
phi1 = sqrt(2/Ts)*cos(2*pi*fc*t);
phi2 = sqrt(2/Ts)*sin(2*pi*fc*t); 
%set snr
SNR_db = -30;
SNR = 10.^(SNR_db/10); %Es/N0 = 20 (db)
N0 = Es./SNR;
%prob of P(si)
P_si = ones(1, M)./M; %when equal prob
%P_si = [4/24 ,13/24, 1/24, 6/24];
%generate S vector (basis:phi1, phi2)
S = [];
for i = 1: M
    si_wave =  sqrt(2*Es/Ts)*cos(2*pi*fc*t- 2*pi*(i-1)/M);
    siXphi1 = si_wave .* phi1; siXphi2 = si_wave .* phi2;
    int_siXphi1 = sample_time*sum(siXphi1);
    int_siXphi2 = sample_time*sum(siXphi2);
    Si = [int_siXphi1, int_siXphi2];
    S = [S ;Si];
end
%generate R vactor (basis:phi1, phi2, all constellation diagram) 
cut_num = 500;
lower_bound = -1.5;
upper_bound = 1.5;
R = [];
add_tern = linspace(lower_bound,upper_bound,cut_num).';
for mul_tern = linspace(lower_bound,upper_bound,cut_num)
    Ri = [ones(cut_num, 1)*mul_tern, add_tern];
    R = [R;Ri];
end
%generate Di
D = zeros(cut_num*cut_num, M);
for i = 1:M
    distance = -abs((R-S(i,:)).^(2));
    D(:, i) = distance(:, 1) + distance(:, 2) + N0*log(P_si(i));
end
D = D.';
%find m_hat
[D_max, m_hat] = max(D);
figure
scatter(R(:,1), R(:,2),[],m_hat.','filled');
title(sprintf('Decision Regions of %d - PSK, SNR = %.2f dB', M, SNR_db))
xlabel('\phi_{1}(t)') 
ylabel('\phi_{2}(t)') 
colormap(summer)
hold on
scatter(S(:,1), S(:,2), [], 'black', 'filled');
text(S(:,1) + 0.1, S(:,2) + 0.1, "S" + string(1:M));
text(-1, -1, "P(S"+ string(1:M)+ ") = {" + string(P_si) + "}");
xlim([-1.5, 1.5])
ylim([-1.5, 1.5])





