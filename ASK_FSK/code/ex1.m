%ex1=================================================
%[1]
%para
A = 1;
f0 = 500;
T = 0.005;
sample_time = 0.0001;
total_time = 0.1;
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
s2_wave = A*sin(2*pi*f0*t);
%s(t)
s_wave = m_wave.*s2_wave;

figure;

%s2(t)
subplot(3, 1, 1);
plot(t, s2_wave);
xlabel('t');
ylabel('s2(t)');

%m(t)
subplot(3, 1, 2);
plot(t, m_wave);
xlabel('t');
ylabel('m(t)');

%s(t)
subplot(3, 1, 3);
plot(t, s_wave);
xlabel('t');
ylabel('m(t)');

%[2]
clear
%para
A = 1;
f1 = 200;
f2 = 600;
T = 0.005;
sample_time = 0.0001;
total_time = 0.1;
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
s1_wave = A*sin(2*pi*f1*t);
%s2(t)
s2_wave = A*sin(2*pi*f2*t);
%s(t)
s_wave = m_wave.*s2_wave;
%s(t)
s_wave = m_wave.*s2_wave + (1-m_wave).*s1_wave;

figure;

%s2(t)
subplot(4, 1, 1);
plot(t, s2_wave);
xlabel('t');
ylabel('s2(t)');

%s2(t)
subplot(4, 1, 2);
plot(t, s1_wave);
xlabel('t');
ylabel('s1(t)');

%m(t)
subplot(4, 1, 3);
plot(t, m_wave);
xlabel('t');
ylabel('m(t)');

%s(t)
subplot(4, 1, 4);
plot(t, s_wave);
xlabel('t');
ylabel('s(t)');

