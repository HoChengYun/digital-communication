%(a)
%A1=20mV
i = 1;
for N0_dB = [-20:-2:-80]
    No_symbols = 100000;
    A = 20*10^(-3);
    T = 1*10^(-3);
    m = round(rand(1, No_symbols)); % symbols 
    s = (2*m-ones(size(m)))*A*T;  %BPSK 
    nois = 10.^(N0_dB/20)*sqrt(T/2);  %noise amplitude 
    n = nois*randn(size(s));  %noise 
    v=s+n; % received signals
    s_hat = (v>0)*2*A*T - A*T;
    ber_20mv(i) = sum((s_hat~=s))/No_symbols;
    ber_20mv_t(i) = Qfunc(sqrt(2*A*A*T/(10.^(N0_dB/10))));
    i = i+1;
end
%A1=30mV
i = 1;
for N0_dB = [-20:-2:-80]
    No_symbols = 100000;
    A = 30*10^(-3);
    T = 1*10^(-3);
    m = round(rand(1, No_symbols)); % symbols 
    s = (2*m-ones(size(m)))*A*T;  %BPSK 
    nois = 10.^(N0_dB/20)*sqrt(T/2);  %noise amplitude 
    n = nois*randn(size(s));  %noise 
    v=s+n; % received signals
    s_hat = (v>0)*2*A*T - A*T;
    ber_30mv(i) = sum((s_hat~=s))/No_symbols;
    ber_30mv_t(i) = Qfunc(sqrt(2*A*A*T/(10.^(N0_dB/10))));
    i = i+1;
end
N0_dB = [-20:-2:-80];
figure
semilogy(N0_dB,ber_20mv,'o',N0_dB,ber_30mv,'o',N0_dB,ber_20mv_t,N0_dB,ber_30mv_t)
ylim([0.000001 1])
set(gca,'XDir','reverse')
title('BER Performance for A=20mV and A=30mV'); 
xlabel('N0 (dB)'); 
ylabel('BER'); 
legend('A=20mV', 'A=30mV', 'A=20mV(theory)', 'A=30mV(theory)', 'Location', 'southwest'); 

%(b)
%T = 1ms
i = 1;
for N0_dB = [-20:-2:-86]
    No_symbols = 100000;
    A = 20*10^(-3);
    T = 1*10^(-3);
    m = round(rand(1, No_symbols)); % symbols 
    s = (2*m-ones(size(m)))*A*T;  %BPSK 
    nois = 10.^(N0_dB/20)*sqrt(T/2);  %noise amplitude 
    n = nois*randn(size(s));  %noise 
    v=s+n; % received signals
    s_hat = (v>0)*2*A*T - A*T;
    ber_1ms(i) = sum((s_hat~=s))/No_symbols;
    ber_1ms_t(i) = Qfunc(sqrt(2*A*A*T/(10.^(N0_dB/10))));
    i = i+1;
end
%T = 0.1ms
i = 1;
for N0_dB = [-20:-2:-86]
    No_symbols = 100000;
    A = 20*10^(-3);
    T = 1*10^(-4);
    m = round(rand(1, No_symbols)); % symbols 
    s = (2*m-ones(size(m)))*A*T;  %BPSK 
    nois = 10.^(N0_dB/20)*sqrt(T/2);  %noise amplitude 
    n = nois*randn(size(s));  %noise 
    v=s+n; % received signals
    s_hat = (v>0)*2*A*T - A*T;
    ber_01ms(i) = sum((s_hat~=s))/No_symbols;
    ber_01ms_t(i) = Qfunc(sqrt(2*A*A*T/(10.^(N0_dB/10))));
    i = i+1;
end
N0_dB = [-20:-2:-86];
figure
semilogy(N0_dB,ber_1ms,'o',N0_dB,ber_01ms,'o',N0_dB,ber_1ms_t,N0_dB,ber_01ms_t)
ylim([0.000001 1])
set(gca,'XDir','reverse')
title('BER Performance for T=1ms and T=0.1ms'); 
xlabel('N0 (dB)'); 
ylabel('BER'); 
legend('T=1ms', 'T=0.1ms','T=1ms(theory)', 'T=0.1ms(theory)', 'Location', 'southwest'); 

%(c)
%A=20mv,T = 1ms
i = 1;
for N0_dB = [-20:-2:-86]
    No_symbols = 100000;
    A = 20*10^(-3);
    T = 1*10^(-3);
    m = round(rand(1, No_symbols)); % symbols 
    s = (2*m-ones(size(m)))*A*T;  %BPSK 
    nois = 10.^(N0_dB/20)*sqrt(T/2);  %noise amplitude 
    n = nois*randn(size(s));  %noise 
    v=s+n; % received signals
    s_hat = (v>0)*2*A*T - A*T;
    ber_ex1a(i) = sum((s_hat~=s))/No_symbols;
    ber_ex1a_t(i) = Qfunc(sqrt(2*A*A*T/(10.^(N0_dB/10))));
    i = i+1;
end
%A=63.2mv,T = 0.1ms
i = 1;
for N0_dB = [-20:-2:-86]
    No_symbols = 100000;
    A = 63.2*10^(-3);
    T = 0.1*10^(-3);
    m = round(rand(1, No_symbols)); % symbols 
    s = (2*m-ones(size(m)))*A*T;  %BPSK 
    nois = 10.^(N0_dB/20)*sqrt(T/2);  %noise amplitude 
    n = nois*randn(size(s));  %noise 
    v=s+n; % received signals
    s_hat = (v>0)*2*A*T - A*T;
    ber_ex1b(i) = sum((s_hat~=s))/No_symbols;
    ber_ex1b_t(i) = Qfunc(sqrt(2*A*A*T/(10.^(N0_dB/10))));
    i = i+1;
end
N0_dB = [-20:-2:-86];
figure
semilogy(N0_dB,ber_ex1a,'o',N0_dB,ber_ex1b,'o',N0_dB,ber_ex1a_t,N0_dB,ber_ex1b_t)
ylim([0.000001 1])
set(gca,'XDir','reverse')
title('BER Performance for ex1 (a) and (b)'); 
xlabel('N0 (dB)'); 
ylabel('BER'); 
legend('(a) A=20mv,T = 1ms', '(b) A=63.2mv,T = 0.1ms','(a) A=20mv,T = 1ms(theory)', '(b) A=63.2mv,T = 0.1ms(theory)', 'Location', 'southwest'); 