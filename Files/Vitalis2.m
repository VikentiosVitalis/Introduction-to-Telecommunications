% Vikentios Vitalis el18803
% fm = 8 + 3 = 11 = 1 + 1 = 2
% am = 3
fm=2000;
am=3;
fs1=20*fm;
Ts1=1/fs1;
duration=1;
Tm=1/fm;
N=floor(duration/Tm);

Amin=0;
Amax=16; % 2 ^ n = 2 ^ 4 = 16 
n=4; % Logw artiohtas ths syxnothtas 
 
% Question a

for i=1:N
    t_samp1(i)=(i-1)*Ts1;
    x_samp1(i)=cos(2*pi*fm*t_samp1(i))*cos(2*pi*(am+2)*fm*t_samp1(i));
end

% Quantizer

delta=(Amax-Amin)/2^n;

partition(1) = Amin+delta/2;
for i=1:1:2^n-1
   partition(i+1) = partition(i)+delta;
end;

x_ind = quantiz(x_samp1,partition);

for i=1:1:N
   if (x_ind(i)+1>2^n)
      x_qnd(i)=partition(2^n)+delta/2;
   else
      x_qnd(i) = partition(x_ind(i)+1)-delta/2;
   end;
end;

k=0;
for i=Amin:0.01:Amax
   k=k+1;
   x_in(k)=i;
   y_ind(k) = quantiz(i,partition);
   if (y_ind(k)+1>2^n)
      y_out(k)=partition(2^n)+delta/2;
   else
      y_out(k)=partition(y_ind(k)+1)-delta/2;
   end;
end;

for i=1:1:k
   if (y_ind(i)<2^n)
      t1=de2bi(y_ind(i),n);
   else
      t1=de2bi(y_ind(i)-1,n);
   end;
   for j=1:1:n
      bin_y_qnd(i,j)=t1(j);
   end;
   
   t2=bin2gray(t1);
   for j=1:1:n
      gray_y_qnd(i,j)=t2(j);
   end;
end;

x_zeros=zeros(1,k)-2.3;
figure(1)
plot(x_in,y_out,'LineWidth',3);
set(gca,'YTick',[]);
text(x_zeros,y_out,num2str(gray_y_qnd),'FontSize',10);
grid;
xlabel(' Input Level ');
title(' Gray Encoded Quantizer I/O ');

figure(2)
plot(t_samp1(1:1:80), x_qnd(1:1:80));
grid;
xlabel('Time (sec)');
ylabel('Quantized Signal Value');
title('Quantizer');

figure(3)
plot(t_samp1(1:1:80),x_samp1(1:1:80),t_samp1(1:1:80),x_qnd(1:1:80));
grid;
xlabel('Time (sec)');
ylabel('Signal Value');
title('Original Signal - Quantized Signal');
legend('original','quantized');
 
for i=1:1:N
   if (x_ind(i)<2^n)
      t1=de2bi(x_ind(i),n);
   else
      t1=de2bi(x_ind(i)-1,n);
   end;
   for j=1:1:n
      bin_x_qnd(i,j)=t1(j);
   end;
   
   t2=bin2gray(t1);
   for j=1:1:n
      gray_x_qnd(i,j)=t2(j);
   end;
end;

figure(4)
stem(t_samp1(1:1:40),x_qnd(1:1:40));
text(t_samp1(1:1:40),x_qnd(1:1:40), num2str(gray_x_qnd(1:40,:)), 'FontSize',5);
grid;
xlabel('Time (sec)');
ylabel('Signal Value');
title('Gray Encoded Signal');

% Question b

q_err=x_samp1-x_qnd;

std_10_samps=std(q_err(1:10));
power_10_samps=sum(x_samp1(1:10).^2)/10;
mse_10_samps=sum(q_err(1:10).^2)/10;
snr_10_samps=power_10_samps/(mse_10_samps/12);
snr_10_samps_db=10*log10(snr_10_samps);

std_20_samps=std(q_err(1:20));
power_20_samps=sum(x_samp1(1:20).^2)/20;
mse_20_samps=sum(q_err(1:20).^2)/20;
snr_20_samps=power_20_samps/(mse_20_samps/12);
snr_20_samps_db=10*log10(snr_20_samps);

disp(' Standard Deviation for 10 samples ');
disp(std_10_samps);

disp(' Standard Deviation for 20 samples ');
disp(std_20_samps);

disp(' SNR for 10 samples ');
disp(snr_10_samps);

disp(' SNR for 20 samples ');
disp(snr_20_samps);

disp(' SNR for 10 samples (db)');
disp(snr_10_samps_db);

disp(' SNR for 20 samples (db)');
disp(snr_20_samps_db);

% Question c

k=0;
for i=1:1:3
   for j=1:1:n
      k=k+1;
      bits(k)=bin_x_qnd(i,j);
   end;
end;

bitrate=n*1000;
Vp=fm/1000;

% Mapping
for i=1:length(bits)
  if (bits(i)==1)
    NRZ_out(i) = Vp;
  else
    NRZ_out(i) = -Vp;
  end;
end;

% Pulse Shaping

i=1;
t=0:0.01:length(bits);
for j=1:length(t)
  if (t(j)<=i)
    y(j)=NRZ_out(i);
  else
    y(j)=NRZ_out(i);
    i=i+1;
  end;
end;

% Plotting
figure(5)
plot(t,y)
axis([0  length(bits) -Vp-4 Vp+4])
grid;
xlabel('Time (sec)');
ylabel('Signal Value');
title('POLAR NRZ');

