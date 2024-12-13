%A = 0.3, 0.5, 0.8 & T = 10^(-3)================
%A = 0.3===========
%para
clear
A = 0.3;
f1 = 1000000;
f2 = 2*1000000;
T = 0.001;
sample_time = 0.0000001;
total_time = 20;
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
s1_wave = A*cos(2*pi*f1*t);
%s2(t)
s2_wave = A*cos(2*pi*f2*t);
%s(t)
s_wave = m_wave.*s2_wave + (1-m_wave).*s1_wave;


%s(t) multi s2_wave-s1_wave = p(t)
p_wave = [];
p_wave = s_wave.*(s2_wave-s1_wave);

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
k_opt = 0;
for N0_i = N0
    sigma_n = sqrt(N0_i/2*A*A*T);
    V =  Q + normrnd(0,sigma_n,1,No_symbols);
    S_hat = V >= k_opt;
    Pe = cat(2,Pe,sum(S_hat~=m)/No_symbols);
    Pe_the = cat(2,Pe_the,Qfunc(sqrt(A.*A*T./(2*N0_i))));
end

Pe_03 = Pe;
Pe_the_03 = Pe_the;

figure;
semilogy(SNR , Pe, "x-", SNR, Pe_the,"o-")
title('FSK BER Performance for A=0.3V & T = 10 (-3) s'); 
xlabel('Eb/N0 (dB)'); 
ylabel('BER'); 
legend('Pe ','(theory)Pe ','Location', 'southwest');


%A = 0.5============
%para
clearvars -except Pe_03 Pe_the_03
A = 0.5;
f1 = 1000000;
f2 = 2*1000000;
T = 0.001;
sample_time = 0.0000001;
total_time = 20;
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
s1_wave = A*cos(2*pi*f1*t);
%s2(t)
s2_wave = A*cos(2*pi*f2*t);
%s(t)
s_wave = m_wave.*s2_wave + (1-m_wave).*s1_wave;


%s(t) multi s2_wave-s1_wave = p(t)
p_wave = [];
p_wave = s_wave.*(s2_wave-s1_wave);

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
k_opt = 0;
for N0_i = N0
    sigma_n = sqrt(N0_i/2*A*A*T);
    V =  Q + normrnd(0,sigma_n,1,No_symbols);
    S_hat = V >= k_opt;
    Pe = cat(2,Pe,sum(S_hat~=m)/No_symbols);
    Pe_the = cat(2,Pe_the,Qfunc(sqrt(A.*A*T./(2*N0_i))));
end

Pe_05 = Pe;
Pe_the_05 = Pe_the;


%A = 0.8============
%para
clearvars -except Pe_03 Pe_the_03 Pe_05 Pe_the_05
A = 0.8;
f1 = 1000000;
f2 = 2*1000000;
T = 0.001;
sample_time = 0.0000001;
total_time = 20;
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
s1_wave = A*cos(2*pi*f1*t);
%s2(t)
s2_wave = A*cos(2*pi*f2*t);
%s(t)
s_wave = m_wave.*s2_wave + (1-m_wave).*s1_wave;


%s(t) multi s2_wave-s1_wave = p(t)
p_wave = [];
p_wave = s_wave.*(s2_wave-s1_wave);

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
k_opt = 0;
for N0_i = N0
    sigma_n = sqrt(N0_i/2*A*A*T);
    V =  Q + normrnd(0,sigma_n,1,No_symbols);
    S_hat = V >= k_opt;
    Pe = cat(2,Pe,sum(S_hat~=m)/No_symbols);
    Pe_the = cat(2,Pe_the,Qfunc(sqrt(A.*A*T./(2*N0_i))));
end

Pe_08 = Pe;
Pe_the_08 = Pe_the;


figure;
semilogy(-N0_db , Pe_08, "x-", -N0_db, Pe_the_08,"o-", -N0_db , Pe_05, "x-", -N0_db, Pe_the_05,"o-",-N0_db , Pe_03, "x-", -N0_db, Pe_the_03,"o-")
title('FSK BER Performance for A=0.8, 0.5, 0.3V & T = 10 (-3) s'); 
xlabel('Eb/N0 (dB)'); 
ylabel('BER'); 
legend('Pe (A = 0.8) ','(theory)Pe (A = 0.8) ','Pe (A = 0.5) ','(theory)Pe (A = 0.5) ','Pe (A = 0.3) ','(theory)Pe (A = 0.3) ','Location', 'southwest');

%[2]
%para
clear
A = 0.3;
f1 = 1000000;
f2 = 2*1000000;
T = 0.0001;
sample_time = 0.0000001;
total_time = 2;
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
s1_wave = A*cos(2*pi*f1*t);
%s2(t)
s2_wave = A*cos(2*pi*f2*t);
%s(t)
s_wave = m_wave.*s2_wave + (1-m_wave).*s1_wave;


%s(t) multi s2_wave-s1_wave = p(t)
p_wave = [];
p_wave = s_wave.*(s2_wave-s1_wave);

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
k_opt = 0;
for N0_i = N0
    sigma_n = sqrt(N0_i/2*A*A*T);
    V =  Q + normrnd(0,sigma_n,1,No_symbols);
    S_hat = V >= k_opt;
    Pe = cat(2,Pe,sum(S_hat~=m)/No_symbols);
    Pe_the = cat(2,Pe_the,Qfunc(sqrt(A.*A*T./(2*N0_i))));
end

figure;
semilogy(SNR , Pe, "x-", SNR, Pe_the,"o-")
title('FSK BER Performance for A=0.3V & T = 10 (-4) s'); 
xlabel('Eb/N0 (dB)'); 
ylabel('BER'); 
legend('Pe ','(theory)Pe ','Location', 'southwest');


%[3]
%A = 0.3====================
%T = 10^(-3)
%para
clear
A = 0.3;
f1 = 1000000;
f2 = 2*1000000;
T = 0.001;
sample_time = 0.0000001;
total_time = 20;
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
s1_wave = A*cos(2*pi*f1*t);
%s2(t)
s2_wave = A*cos(2*pi*f2*t);
%s(t)
s_wave = m_wave.*s2_wave + (1-m_wave).*s1_wave;


%s(t) multi s2_wave-s1_wave = p(t)
p_wave = [];
p_wave = s_wave.*(s2_wave-s1_wave);

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
k_opt = 0;
for N0_i = N0
    sigma_n = sqrt(N0_i/2*A*A*T);
    V =  Q + normrnd(0,sigma_n,1,No_symbols);
    S_hat = V >= k_opt;
    Pe = cat(2,Pe,sum(S_hat~=m)/No_symbols);
    Pe_the = cat(2,Pe_the,Qfunc(sqrt(A.*A*T./(2*N0_i))));
end

Pe_103 = Pe;
Pe_the_103 = Pe_the;

%T = 10^(-4)
%para
clearvars -except Pe_103 Pe_the_103
A = 0.3;
f1 = 1000000;
f2 = 2*1000000;
T = 0.0001;
sample_time = 0.0000001;
total_time = 2;
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
s1_wave = A*cos(2*pi*f1*t);
%s2(t)
s2_wave = A*cos(2*pi*f2*t);
%s(t)
s_wave = m_wave.*s2_wave + (1-m_wave).*s1_wave;


%s(t) multi s2_wave-s1_wave = p(t)
p_wave = [];
p_wave = s_wave.*(s2_wave-s1_wave);

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
k_opt = 0;
for N0_i = N0
    sigma_n = sqrt(N0_i/2*A*A*T);
    V =  Q + normrnd(0,sigma_n,1,No_symbols);
    S_hat = V >= k_opt;
    Pe = cat(2,Pe,sum(S_hat~=m)/No_symbols);
    Pe_the = cat(2,Pe_the,Qfunc(sqrt(A.*A*T./(2*N0_i))));
end

Pe_104 = Pe;
Pe_the_104 = Pe_the;

figure;
semilogy(-N0_db , Pe_103, "x-", -N0_db, Pe_104,"o-")
title('Comparisons of BERs of FSK A = 0.3'); 
xlabel('-N0 (dB)'); 
ylabel('BER'); 
legend('T = 10 (-3) ','T = 10 (-4) ','Location', 'southwest');




