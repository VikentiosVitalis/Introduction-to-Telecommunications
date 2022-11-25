% Vikentios Vitalis el18803
% fm = 8 + 3 = 11 = 1 + 1 = 2
% am = 3
fm=2000;		
am = 3;
fs1=20*fm;
fs2=100*fm;
fs3=5*fm;
Ts1=1/fs1;
Ts2=1/fs2;
Ts3=1/fs3;
dur=1;
Tm=1/fm;
N=dur/Tm;

% Question a(i) - Sample with fs1 = 20fm

for i=1:N
   t_samp1(i)=(i-1)*Ts1;
   x_samp1(i)=cos(2*pi*fm*t_samp1(i))*cos(2*pi*(am+2)*fm*t_samp1(i));
end

% Question a(ii) - Sample with fs2 = 100fm

for i=1:N
   t_samp2(i)=(i-1)*Ts2;
   x_samp2(i)=cos(2*pi*fm*t_samp2(i))*cos(2*pi*(am+2)*fm*t_samp2(i));
end

figure(1)
stem(t_samp1(1:41),x_samp1(1:41));
grid;
xlabel('Time Axis');
ylabel('Signal');
title('Sampling Frequency fs1=20fm');

figure(2)
stem(t_samp2(1:201),x_samp2(1:201));
grid;
xlabel('Time Axis');
ylabel('Signal');
title('Sampling Frequency fs2=100fm');

figure(3)
stem(t_samp2(1:201),x_samp2(1:201),'g');
hold on;
stem(t_samp1(1:41),x_samp1(1:41),'k');
grid;
xlabel(' Time axis(sec) ');
ylabel(' Signal ');
title('Common Graph Of fs1 and fs2');
legend('fs1=20*fm','fs2=100*fm');
hold off

% Question b with Fs3=5fm

for i=1:N
   t_samp3(i)=(i-1)*Ts3;
   x_samp3(i)=cos(2*pi*fm*t_samp3(i))*cos(2*pi*(am+2)*fm*t_samp3(i));
end

hf_fs=500*fm;
hf_Ts=1/hf_fs;
hf_N=floor(Tm/hf_Ts);
for i=1:1:hf_N
    t_samp_hf(i)=(i-1)*hf_Ts;
    x_samp_hf(i)=cos(2*pi*fm*t_samp_hf(i))*cos(2*pi*(am+2)*fm*t_samp_hf(i));
end
x_fft_meas=abs(fft(x_samp_hf));
x_fft_meas_db=20*log10(abs(fft(x_samp_hf)));

% Graphs

figure(4)
stem(t_samp3(1:N/100),x_samp3(1:N/100));
grid;
xlabel('Time Axis');
ylabel('Signal');
title('Sampling Frequency fs3=5fm');

figure(5)
i=1:1:hf_N/2;
plot((i-1)*hf_fs/hf_N/1000,x_fft_meas_db(i),'LineWidth',2);
grid;
xlabel('Frequency (KHz)');
ylabel('FFT Amplitude Value (db)');
title('Fourier Analysis of Input Signal');