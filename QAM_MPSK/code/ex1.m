%para
clear
fc = 3*10^(6);
Ts = 10^(-6);
Es = 10^(-6);
sample_frequency = 10^(10);
sample_time = 1/sample_frequency;
No_symbols = 8;
total_time = No_symbols*Ts;
t = [0:sample_time:total_time-sample_time];
phi1 = sqrt(2/Ts)*cos(2*pi*fc*t);
phi2 = sqrt(2/Ts)*sin(2*pi*fc*t);

s = [sqrt(Es/2)*phi1+sqrt(Es/2)*phi2; -sqrt(Es/2)*phi1+sqrt(Es/2)*phi2;...
    -sqrt(Es/2)*phi1-sqrt(Es/2)*phi2; sqrt(Es/2)*phi1-sqrt(Es/2)*phi2];
%generate m
%m = randi([1,4],1,No_symbols);
m = [1 2 3 4 1 2 3 4];
s_wave = [];
tep_point = 1;
Ts_point = round(Ts/sample_time); 
for i = m
    si = s(i,:);
    wave_part = si(tep_point:tep_point+Ts_point-1);
    s_wave = cat(2, s_wave, wave_part);
    tep_point = tep_point+Ts_point;
end

plot(t,s_wave)

