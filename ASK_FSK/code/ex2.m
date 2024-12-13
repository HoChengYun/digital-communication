%ex2=====================================================
%(1.)------------------------------------------
%[1]
%T = 10^(-3)
clear
T = 0.001;
A = 0.2;
N0_db = [-15:-3:-60];
N0 = 10.^(N0_db/10);
Eb = A*A*(T/2);
Pe = [];
for i = N0
    Pe = cat(2,Pe,Qfunc(sqrt(A.*A*T/(4*i))));
end
Pe1 = Pe;
N01 = N0;
Eb1 = Eb;
clearvars -except Pe1 N01 Eb1

%T = 10^(-4)
T = 0.0001;
A = 0.2;
N0_db = [-25:-3:-70];
N0 = 10.^(N0_db/10);
Eb = A*A*(T/2);
Pe = [];
for i = N0
    Pe = cat(2,Pe,Qfunc(sqrt(A.*A*T/(4*i))));
end
Pe2 = Pe;
N02 = N0;
Eb2 = Eb;

SNR1 = 10*log10(Eb1./N01);
SNR2 = 10*log10(Eb2./N02);
figure;
semilogy(SNR1 , Pe1, "x-", SNR2, Pe2,"o-")
title('BER Performance for T=1ms and T=0.1ms'); 
xlabel('Eb/N0 (dB)'); 
ylabel('BER'); 
legend('Pe for Tb1 = 10 (-3)','Pe for Tb2 = 10 (-4)','Location', 'southwest');

%[2]
%para
clear
A = 0.2;
f0 = 10000;
T = 0.001;
sample_time = 0.00001;
total_time = 10;
t = [0:sample_time:total_time-sample_time];
No_symbols = total_time/T;
m = round(rand(1, No_symbols));
%generate m(t), s(t), s1(t), s2(t)
m_wave = [];
s_wave = [];
s1_wave = [];
s2_wave = [];
s_wave = [];
%m(t)
for i = m
    m_wave_part = ones(1,T/sample_time)*i;
    m_wave = cat(2, m_wave, m_wave_part);
end
%s1(t)
s1_wave = [];
%s2(t)
s2_wave = A*cos(2*pi*f0*t);
%s(t)
s_wave = m_wave.*s2_wave;

%s(t) multi Acos(wct) = p(t)
p_wave = [];
p_wave = s_wave.*A.*cos(2*pi*f0*t);

%integral [0, T] = Q
P = reshape(p_wave,[T/sample_time,No_symbols]);
Q = sample_time.*sum(P);

% plus noise = V & S_hat estimate
V = [];
Pe = [];
Pe_the = [];
N0_db = [-15:-3:-60];
N0 = 10.^(N0_db/10);
Eb = A*A*(T/2);
SNR = 10*log10(Eb./N0);
k_opt = A*A*T/4;
for N0_i = N0
    sigma_n = sqrt(N0_i/4*A*A*T);
    V =  Q + normrnd(0,sigma_n,1,No_symbols);
    S_hat = V >= k_opt;
    Pe = cat(2,Pe,sum(S_hat~=m)/No_symbols);
    Pe_the = cat(2,Pe_the,Qfunc(sqrt(A.*A*T./(4*N0_i))));
end

figure;
semilogy(SNR , Pe, "x-", SNR, Pe_the,"o-")
title('BER Performance for T=1ms'); 
xlabel('Eb/N0 (dB)'); 
ylabel('BER'); 
legend('Pe for Tb1 = 10 (-3)','(theory)Pe for Tb1 = 10^ (-3)','Location', 'southwest');

%(2.)------------------------------------------
% A = 0.8
%para
clear
A = 0.8;
f0 = 10000;
T = 0.0001;
sample_time = 0.00001;
total_time = 10;
t = [0:sample_time:total_time-sample_time];
No_symbols = total_time/T;
m = round(rand(1, No_symbols));
%generate m(t), s(t), s1(t), s2(t)
m_wave = [];
s_wave = [];
s1_wave = [];
s2_wave = [];
s_wave = [];
%m(t)
for i = m
    m_wave_part = ones(1,round(T/sample_time))*i;
    m_wave = cat(2, m_wave, m_wave_part);
end
%s1(t)
s1_wave = [];
%s2(t)
s2_wave = A*cos(2*pi*f0*t);
%s(t)
s_wave = m_wave.*s2_wave;

%s(t) multi Acos(wct) = p(t)
p_wave = [];
p_wave = s_wave.*A.*cos(2*pi*f0*t);

%integral [0, T] = Q
P = reshape(p_wave,[round(T/sample_time),No_symbols]);
Q = sample_time.*sum(P);

% plus noise = V & S_hat estimate
V = [];
Pe = [];
Pe_the = [];
N0_db = [-15:-3:-60];
N0 = 10.^(N0_db/10);
Eb = A*A*(T/2);
SNR = 10*log10(Eb./N0);
k_opt = A*A*T/4;
for N0_i = N0
    sigma_n = sqrt(N0_i/4*A*A*T);
    V =  Q + normrnd(0,sigma_n,1,No_symbols);
    S_hat = V >= k_opt;
    Pe = cat(2,Pe,sum(S_hat~=m)/No_symbols);
    Pe_the = cat(2,Pe_the,Qfunc(sqrt(A.*A*T./(4*N0_i))));
end

figure;
semilogy(SNR , Pe, "x-", SNR, Pe_the,"o-")
title('BER Performance for A = 0.8'); 
xlabel('Eb/N0 (dB)'); 
ylabel('BER'); 
legend('Pe for A1 = 0.8','(theory)Pe for A1 = 0.8','Location', 'southwest');
Pe_08 = Pe;
Pe_the_08 = Pe_the;


% A = 0.5
%para
clearvars -except Pe_08 Pe_the_08
A = 0.5;
f0 = 10000;
T = 0.0001;
sample_time = 0.00001;
total_time = 10;
t = [0:sample_time:total_time-sample_time];
No_symbols = total_time/T;
m = round(rand(1, No_symbols));
%generate m(t), s(t), s1(t), s2(t)
m_wave = [];
s_wave = [];
s1_wave = [];
s2_wave = [];
s_wave = [];
%m(t)
for i = m
    m_wave_part = ones(1,round(T/sample_time))*i;
    m_wave = cat(2, m_wave, m_wave_part);
end
%s1(t)
s1_wave = [];
%s2(t)
s2_wave = A*cos(2*pi*f0*t);
%s(t)
s_wave = m_wave.*s2_wave;

%s(t) multi Acos(wct) = p(t)
p_wave = [];
p_wave = s_wave.*A.*cos(2*pi*f0*t);

%integral [0, T] = Q
P = reshape(p_wave,[round(T/sample_time),No_symbols]);
Q = sample_time.*sum(P);

% plus noise = V & S_hat estimate
V = [];
Pe = [];
Pe_the = [];
N0_db = [-15:-3:-60];
N0 = 10.^(N0_db/10);
Eb = A*A*(T/2);
SNR = 10*log10(Eb./N0);
k_opt = A*A*T/4;
for N0_i = N0
    sigma_n = sqrt(N0_i/4*A*A*T);
    V =  Q + normrnd(0,sigma_n,1,No_symbols);
    S_hat = V >= k_opt;
    Pe = cat(2,Pe,sum(S_hat~=m)/No_symbols);
    Pe_the = cat(2,Pe_the,Qfunc(sqrt(A.*A*T./(4*N0_i))));
end

figure;
semilogy(SNR , Pe, "x-", SNR, Pe_the,"o-")
title('BER Performance for A = 0.5'); 
xlabel('Eb/N0 (dB)'); 
ylabel('BER'); 
legend('Pe for A1 = 0.5','(theory)Pe for A1 = 0.5','Location', 'southwest');
Pe_05 = Pe;
Pe_the_05 = Pe_the;


%A = 0.8 compare A = 0.5
figure;
semilogy(-N0_db , Pe_05, "x-", -N0_db, Pe_08,"o-")
title('BER Performance for A1 = 0.5 & A2 = 0.8'); 
xlabel('-N0 (dB)'); 
ylabel('BER'); 
legend('Pe for A1 = 0.5','Pe for A2 = 0.8','Location', 'southwest');

