%para
clear
fc = 3*10^(6);
Ts = 10^(-6);
a = 10^(-3);
sample_frequency = 10^(8);
sample_time = 1/sample_frequency;
No_symbols = 200;
total_time = No_symbols*Ts;
t = [0:sample_time:total_time-sample_time];
phi1 = sqrt(2/Ts)*cos(2*pi*fc*t);
phi2 = sqrt(2/Ts)*sin(2*pi*fc*t); 

trial = 0;
all_s_wave = [];
for trial = 1:300
    select_1_3 = [1 3 -1 -3];
    m = reshape(randsample(select_1_3, No_symbols*2, true), No_symbols, 2);
    s_wave = [];
    tep_point = 1;
    Ts_point = round(Ts/sample_time); 
    for i = 1: No_symbols
        si = a*m(i,1)*phi1 + a*m(i,2)*phi2;
        wave_part = si(tep_point:tep_point+Ts_point-1);
        s_wave = cat(2, s_wave, wave_part);
        tep_point = tep_point+Ts_point;
    end
    all_s_wave = [all_s_wave; s_wave];
end
all_s_wave = all_s_wave'; %columns: different trials

%para set
Es =  10*a*a; %average energy
t_one = -Ts/2:sample_time:Ts/2-sample_time;
N = length(t_one);
N_zero = 50*N;
df = sample_frequency/(N+N_zero);
f_one = -sample_frequency/2:df:sample_frequency/2-df;
%make time basis func
phi1_one = sqrt(2/Ts)*cos(2*pi*fc*t_one);
phi2_one = sqrt(2/Ts)*sin(2*pi*fc*t_one);
%make frequency basis func
PHI1_one = fftshift(fft([phi1_one zeros(1,N_zero)]))/sample_frequency;
PHI2_one = fftshift(fft([phi2_one zeros(1,N_zero)]))/sample_frequency;

%theory psd
theory_psd = Es/(2*Ts) * (PHI1_one.*conj(PHI1_one) + PHI2_one.*conj(PHI2_one));



%experiment psd
N = 5*length(t);
df = sample_frequency/N;
f = -sample_frequency/2:df:sample_frequency/2-df;
all_S_WAVE = fftshift(fft(all_s_wave, N))/(sample_frequency);
experiment_psd = 1/(No_symbols*Ts) * all_S_WAVE .* conj(all_S_WAVE);
average_psd = sum(experiment_psd,2)/trial;

figure
plot(f,experiment_psd)
xlim([-1.5*10^(7) 1.5*10^(7)])
title(sprintf('%d overlapped spectra for Ts = %d μs', trial, Ts)); 
xlabel('freq(Hz)'); 
ylabel('PSD(W/Hz)'); 

figure
plot(f, average_psd, '--', 'LineWidth', 1)
hold on;
plot(f_one,theory_psd, 'LineWidth', 1)
xlim([-1.5*10^(7) 1.5*10^(7)])
title('psd for Ts = ', num2str(Ts)); 
xlabel('freq(Hz)'); 
ylabel('PSD(W/Hz)'); 
legend('experiment','theory','Location', 'northeast');



%para
clear
fc = 3*10^(6);
Ts = 4*10^(-6);
a = 10^(-3);
sample_frequency = 10^(8);
sample_time = 1/sample_frequency;
No_symbols = 200;
total_time = No_symbols*Ts;
t = [0:sample_time:total_time-sample_time];
phi1 = sqrt(2/Ts)*cos(2*pi*fc*t);
phi2 = sqrt(2/Ts)*sin(2*pi*fc*t); 

trial = 0;
all_s_wave = [];
for trial = 1:300
    select_1_3 = [1 3 -1 -3];
    m = reshape(randsample(select_1_3, No_symbols*2, true), No_symbols, 2);
    s_wave = [];
    tep_point = 1;
    Ts_point = round(Ts/sample_time); 
    for i = 1: No_symbols
        si = a*m(i,1)*phi1 + a*m(i,2)*phi2;
        wave_part = si(tep_point:tep_point+Ts_point-1);
        s_wave = cat(2, s_wave, wave_part);
        tep_point = tep_point+Ts_point;
    end
    all_s_wave = [all_s_wave; s_wave];
end
all_s_wave = all_s_wave'; %columns: different trials

%para set
Es =  10*a*a; %average energy
t_one = -Ts/2:sample_time:Ts/2-sample_time;
N = length(t_one);
N_zero = 50*N;
df = sample_frequency/(N+N_zero);
f_one = -sample_frequency/2:df:sample_frequency/2-df;
%make time basis func
phi1_one = sqrt(2/Ts)*cos(2*pi*fc*t_one);
phi2_one = sqrt(2/Ts)*sin(2*pi*fc*t_one);
%make frequency basis func
PHI1_one = fftshift(fft([phi1_one zeros(1,N_zero)]))/sample_frequency;
PHI2_one = fftshift(fft([phi2_one zeros(1,N_zero)]))/sample_frequency;

%theory psd
theory_psd = Es/(2*Ts) * (PHI1_one.*conj(PHI1_one) + PHI2_one.*conj(PHI2_one));



%experiment psd
N = 5*length(t);
df = sample_frequency/N;
f = -sample_frequency/2:df:sample_frequency/2-df;
all_S_WAVE = fftshift(fft(all_s_wave, N))/(sample_frequency);
experiment_psd = 1/(No_symbols*Ts) * all_S_WAVE .* conj(all_S_WAVE);
average_psd = sum(experiment_psd,2)/trial;

figure
plot(f, average_psd, '--', 'LineWidth', 1)
hold on;
plot(f_one,theory_psd, 'LineWidth', 1)
xlim([-5*10^(6) 5*10^(6)])
title('psd for Ts = ', num2str(Ts)); 
xlabel('freq(Hz)'); 
ylabel('PSD(W/Hz)'); 
legend('experiment','theory','Location', 'north');




