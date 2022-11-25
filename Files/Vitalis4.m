% Vikentios Vitalis el18803
% 8 + 3 = 11 = 1 + 1 = 2
% shannon_even.txt logw artiothtas
clc;
clear;

% Question a 
   nbits = 46;
   randomseq = (1+sign(randn(1,nbits)))/2;
   amplitude = 2;
   samps_per_bit = 16;
   Tbit = 0.5;
   fc = 1/Tbit;


%	QPSK

   v=0;
   for i=1:1:nbits/2
      for j=1:1:samps_per_bit
         v=v+1;
         timepoint(v)=(Tbit/samps_per_bit)*(v-1);
         if (randomseq(2*(i-1)+1:2*(i-1)+2)==00)
            x_sign(v)=amplitude*cos(2*pi*fc*timepoint(v)+pi/4);
            d_sign(v)=0;
         elseif (randomseq(2*(i-1)+1:2*(i-1)+2)==01)
            x_sign(v)=amplitude*cos(2*pi*fc*timepoint(v)+pi/4+3*pi/2);
            d_sign(v)=1;
         elseif (randomseq(2*(i-1)+1:2*(i-1)+2)==10)
            x_sign(v)=amplitude*cos(2*pi*fc*timepoint(v)+pi/4+pi/2);
            d_sign(v)=2;
         else
            x_sign(v)=amplitude*cos(2*pi*fc*timepoint(v)+pi/4+pi);
            d_sign(v)=3;
         end;
      end;
   end;
   
   L=4;
   % Mapping vector for M - PSK Gray encoding
   k=log2(L); 			% Number of bits per point

   ph1 = [pi/4]; 
   theta1 = [ph1; -ph1; pi-ph1; -pi+ph1]; 
   map = exp(1j*theta1); 
   if(k>2) 
      for j=3:k 
         theta1=theta1/2; 
         map=exp(1j*theta1); 
         map=[map; -conj(map)]; 
         theta1=real(log(map)/1j); 
      end;
   end;
   xstar_nonoise=(amplitude^2/2)^0.5*real(map);
   ystar_nonoise=(amplitude^2/2)^0.5*imag(map);
  
   figure(1)
   plot(xstar_nonoise,ystar_nonoise,'x');
   for i=1:1:L
      text(xstar_nonoise(i)-0.5,ystar_nonoise(i)-0.5,num2str(de2bi(i-1,log2(L),'left-msb')), 'FontSize', 5);
   end;
   grid;
   axis([-max(xstar_nonoise)-1 max(xstar_nonoise)+1 -max(xstar_nonoise)-1 max(xstar_nonoise)+1]);
   xlabel('(Eb/Tb)^0^.^5 (I)');
   ylabel('(Eb/Tb)^0^.^5 (Q)');
   title('Constellation Map for Q-PSK');

   figure(2)
   plot(timepoint,x_sign,'LineWidth',1.2);
   hold on;
   plot(timepoint,d_sign,'LineWidth',1.2);
   hold off;
   for i=1:samps_per_bit:length(d_sign)
      text(timepoint(i)+0.1,d_sign(i)+0.5,num2str(de2bi(d_sign(i),log2(L),'left-msb')), 'FontSize', 5);
   end;
   grid;
   axis([0 max(timepoint)+1 -1-max(x_sign) 1+max(x_sign)]);
   xlabel('Time (s)');
   ylabel('Signal Value');
   title('Q-PSK Signal');

% Question b

   SNR=5;

   xsymb=bi2de(reshape(randomseq,log2(L),length(randomseq)/log2(L)).','left-msb'); 
   x_map=[]; 
   for k=1:length(xsymb) 
      x_val(k)=real(map(xsymb(k)+1));
      y_val(k)=imag(map(xsymb(k)+1));
      x_map(k)=x_val(k)+exp(pi/2i)*y_val(k); 
   end;
   x_map=(amplitude^2/2)^0.5*x_map;

   x_sign_energy=sum(abs(x_map).^2);
   
   noise_sign_energy=x_sign_energy/10^(SNR/10);
   noise_length=length(x_map);
   sigma_n=(noise_sign_energy/noise_length)^0.5;
   noise_sigr=normrnd(0,sigma_n,[1,length(x_map)]);
   noise_sigi=normrnd(0,sigma_n,[1,length(x_map)]);
   noise_sign=noise_sigr+noise_sigi*exp(pi/2i);
   noise_sign_energy_calc=sum(abs(noise_sign).^2);
   
   SNRcalculation=10*log10(x_sign_energy/noise_sign_energy_calc);
   SNRresult=SNR-SNRcalculation;
   noise_sign=10^(-SNRresult/10)*noise_sign;
   
   x_sign_n=x_map+noise_sign;

   xstar=real(x_sign_n);
   ystar=imag(x_sign_n);
   
   figure(3)
   plot(xstar,ystar,'x');
   hold on
   plot(xstar_nonoise,ystar_nonoise,'o');
   hold off
   grid;
   axis([-max(xstar)-1 max(xstar)+1 -max(xstar)-1 max(xstar)+1]);
   xlabel('(Eb/Tb)^0.5 (I)');
   ylabel('(Eb/Tb)^0.5 (Q)');
   title('Constellation Map for QPSK (Eb/No=5 db)');
   legend('x=Stars for signal with noise','o=Stars for noiseless signal');



   SNR=15;

   xsymb=bi2de(reshape(randomseq,log2(L),length(randomseq)/log2(L)).','left-msb'); 
   x_map=[]; 
   for k=1:length(xsymb) 
      x_val(k)=real(map(xsymb(k)+1));
      y_val(k)=imag(map(xsymb(k)+1));
      x_map(k)=x_val(k)+exp(pi/2i)*y_val(k); 
   end;
   x_map=(amplitude^2/2)^0.5*x_map;

   x_sign_energy=sum(abs(x_map).^2);
   
   noise_sign_energy=x_sign_energy/10^(SNR/10);
   noise_length=length(x_map);
   sigma_n=(noise_sign_energy/noise_length)^0.5;
   noise_sigr=normrnd(0,sigma_n,[1,length(x_map)]);
   noise_sigi=normrnd(0,sigma_n,[1,length(x_map)]);
   noise_sign=noise_sigr+noise_sigi*exp(pi/2i);
   noise_sign_energy_calc=sum(abs(noise_sign).^2);
   
   SNRcalculation=10*log10(x_sign_energy/noise_sign_energy_calc);
   SNRresult=SNR-SNRcalculation;
   noise_sign=10^(-SNRresult/10)*noise_sign;
   
   x_sign_n=x_map+noise_sign;

   xstar=real(x_sign_n);
   ystar=imag(x_sign_n);
   
   figure(4)
   plot(xstar,ystar,'x');
   hold on
   plot(xstar_nonoise,ystar_nonoise,'o');
   hold off
   grid;
   axis([-max(xstar)-1 max(xstar)+1 -max(xstar)-1 max(xstar)+1]);
   xlabel('(Eb/Tb)^0.5 (I)');
   ylabel('(Eb/Tb)^0.5 (Q)');
   title('Constellation Map for QPSK (Eb/No=15 db)');
   legend('x=Stars for signal with noise','o=Stars for noiseless signal');
   
% Question c - BER vs SNR

   format long;
   
   SNR_start=0;
   SNR_stop=15;
   SNR_step=1;
   SNR_iter=floor((SNR_stop-SNR_start)/SNR_step)+1;

   nbits=500000;
   Tbit=0.5;

   randomseq=(1+sign(randn(nbits,1)))/2;

   amplitude=9;
   samps_per_bit=1;
   
   xsymb=bi2de(reshape(randomseq,log2(L),length(randomseq)/log2(L)).','left-msb'); 
   x_map=[]; 
   for k=1:length(xsymb) 
      x_val(k)=real(map(xsymb(k)+1));
      y_val(k)=imag(map(xsymb(k)+1));
      x_map(k)=x_val(k)+exp(pi/2i)*y_val(k); 
   end;
   x_map=(amplitude^2/2)^0.5*x_map;

   x_sign_energy=sum(abs(x_map).^2);
   
   for s=1:1:SNR_iter

      SNR_st(s)=SNR_start+(s-1)*SNR_step;
      SNR=SNR_st(s);
   
      noise_sign_energy=x_sign_energy/10^(SNR/10);
      noise_length=length(x_map);
      sigma_n=(noise_sign_energy/noise_length)^0.5;
      noise_sigr=normrnd(0,sigma_n,[1,length(x_map)]);
      noise_sigi=normrnd(0,sigma_n,[1,length(x_map)]);
      noise_sign=noise_sigr+noise_sigi*exp(pi/2i);
      noise_sign_energy_calc=sum(abs(noise_sign).^2);
   
      SNRcalculation=10*log10(x_sign_energy/noise_sign_energy_calc);
      SNRresult=SNR-SNRcalculation;
      noise_sign=10^(-SNRresult/10)*noise_sign;
   
      x_sign_n=x_map+noise_sign;

      for i=1:1:nbits/2
         dist0=abs(x_sign_n(i)-amplitude*(1+exp(pi/2i)));
         dist1=abs(x_sign_n(i)-amplitude*(1-exp(pi/2i)));
         dist2=abs(x_sign_n(i)-amplitude*(-1+exp(pi/2i)));
         dist3=abs(x_sign_n(i)-amplitude*(-1-exp(pi/2i)));
         if (min(min(min(dist0,dist1),dist2),dist3)==dist0)
            randomseq_rx(1+2*(i-1))=0;
            randomseq_rx(2+2*(i-1))=0;
         elseif (min(min(min(dist0,dist1),dist2),dist3)==dist1)
            randomseq_rx(1+2*(i-1))=0;
            randomseq_rx(2+2*(i-1))=1;
         elseif (min(min(min(dist0,dist1),dist2),dist3)==dist2)
            randomseq_rx(1+2*(i-1))=1;
            randomseq_rx(2+2*(i-1))=0;
         else
            randomseq_rx(1+2*(i-1))=1;
            randomseq_rx(2+2*(i-1))=1;
         end;
      end;
   
      BER(s)=0;
      for i=1:1:nbits
         if (randomseq(i)~=randomseq_rx(i))
            BER(s)=BER(s)+1;
         end;
      end;
      BER_st(s)=BER(s)/nbits;
            
   end;
   
   figure(5)
   semilogy(SNR_st,BER_st,'LineWidth',2);
   grid;
   xlabel('SNR (db)');
   ylabel('BER');
   title('BER vs SNR for Q-PSK');


% Question d

%(i) 
   filename='shannon_even.txt';
   fid=fopen(filename);
   st_struct_str=fscanf(fid,'%c');
   ld=fclose(fid);

   disp('---------------------------------------');
   disp(' Send ASCII ');
   disp('---------------------------------------');
   disp(st_struct_str)
   
   st_struct=double(st_struct_str);

   sym_bits=8;
   for i=1:1:length(st_struct)
      binfile(i,:)=de2bi(st_struct(i),sym_bits,'left-msb');
   end;
   
   v=0;
   for i=1:1:length(st_struct)
      for j=1:1:sym_bits
         v=v+1;
         bin_stream(v)=binfile(i,j);
      end;
   end;
   
   figure(6)
   i=1:1:length(st_struct);
   semilogy(i,st_struct,'LineWidth',1.5);
   grid;
   xlabel('Sample No');
   ylabel('Value');
   title('Original ASCII Text Signal');
   

%(ii)
   Amax=max(st_struct);
   Amin=min(st_struct);
   n=8;
   N=length(st_struct);

   delta=(Amax-Amin)/2^n;

   partition(1) = Amin+delta/2;
   for i=1:1:2^n-1
      partition(i+1)=partition(i)+delta;
   end;

   x_ind = quantiz(st_struct,partition);

   for i=1:1:N
      if (x_ind(i)+1>2^n)
         x_qnd(i)=partition(2^n)+delta/2;
      else
         x_qnd(i)=partition(x_ind(i)+1);
      end;
   end;

   x_bin=de2bi(floor(x_qnd),sym_bits,'left-msb');

   v=0;
   for i=1:1:length(st_struct)
      for j=1:1:sym_bits
         v=v+1;
         x_bin_str(v)=x_bin(i,j);
      end;
   end;

   figure(7)
   i=1:1:length(x_qnd);
   semilogy(i,x_qnd,'LineWidth',1.5);
   grid;
   xlabel('Sample No');
   ylabel('Value');
   title('Quantized ASCII Text Signal');
   
%(iii-iv)

   SNR=5;
   amplitude=1;

   xsymb=bi2de(reshape(x_bin_str,log2(L),length(x_bin_str)/log2(L)).','left-msb'); 
   x_map=[]; 
   for k=1:length(xsymb) 
      x_val(k)=real(map(xsymb(k)+1));
      y_val(k)=imag(map(xsymb(k)+1));
      x_map(k)=x_val(k)+exp(pi/2i)*y_val(k); 
   end;
   x_map=(amplitude^2/2)^0.5*x_map;

   x_sign_energy=sum(abs(x_map).^2);
   
   noise_sign_energy=x_sign_energy/10^(SNR/10);
   noise_length=length(x_map);
   sigma_n=(noise_sign_energy/noise_length)^0.5;
   noise_sigr=normrnd(0,sigma_n,[1,length(x_map)]);
   noise_sigi=normrnd(0,sigma_n,[1,length(x_map)]);
   noise_sign=noise_sigr+noise_sigi*exp(pi/2i);
   noise_sign_energy_calc=sum(abs(noise_sign).^2);
   
   SNRcalculation=10*log10(x_sign_energy/noise_sign_energy_calc);
   SNRresult=SNR-SNRcalculation;
   noise_sign=10^(-SNRresult/10)*noise_sign;
   
   x_sign_n_5=x_map+noise_sign;

   xstar=real(x_sign_n_5);
   ystar=imag(x_sign_n_5);
   xstar_nonoise=real(map);
   ystar_nonoise=imag(map);
   
   figure(8)
   plot(xstar,ystar,'x');
   hold on
   plot(xstar_nonoise,ystar_nonoise,'o');
   hold off
   grid;
   axis([-max(xstar)-1 max(xstar)+1 -max(xstar)-1 max(xstar)+1]);
   xlabel('(Eb/Tb)^0.5 (I)');
   ylabel('(Eb/Tb)^0.5 (Q)');
   title('Constellation Map for Q-PSK (Eb/No=5 db)');
   legend('x=Stars for signal with noise','o=Stars for noiseless signal');

   SNR = 15;
   
   noise_sign_energy=x_sign_energy/10^(SNR/10);
   noise_length=length(x_map);
   sigma_n=(noise_sign_energy/noise_length)^0.5;
   noise_sigr=normrnd(0,sigma_n,[1,length(x_map)]);
   noise_sigi=normrnd(0,sigma_n,[1,length(x_map)]);
   noise_sign=noise_sigr+noise_sigi*exp(pi/2i);
   noise_sign_energy_calc=sum(abs(noise_sign).^2);
   
   SNRcalculation=10*log10(x_sign_energy/noise_sign_energy_calc);
   SNRresult=SNR-SNRcalculation;
   noise_sign=10^(-SNRresult/10)*noise_sign;
   
   x_sign_n_15=x_map+noise_sign;

   xstar=real(x_sign_n_15);
   ystar=imag(x_sign_n_15);
   xstar_nonoise=real(map);
   ystar_nonoise=imag(map);
   
   figure(9)
   plot(xstar,ystar,'x');
   hold on
   plot(xstar_nonoise,ystar_nonoise,'o');
   hold off
   grid;
   axis([-max(xstar)-1 max(xstar)+1 -max(xstar)-1 max(xstar)+1]);
   xlabel('(Eb/Tb)^0.5 (I)');
   ylabel('(Eb/Tb)^0.5 (Q)');
   title('Constellation Map for Q-PSK (Eb/No=15 db)');
   legend('x=Stars for signal with noise','o=Stars for noiseless signal');

% 	(v) - (vi)

   SNR = 5;

   for i=1:1:length(x_sign_n_5)
      dist0=abs(x_sign_n_5(i)-amplitude*(1+exp(pi/2i)));
      dist1=abs(x_sign_n_5(i)-amplitude*(1-exp(pi/2i)));
      dist2=abs(x_sign_n_5(i)-amplitude*(-1+exp(pi/2i)));
      dist3=abs(x_sign_n_5(i)-amplitude*(-1-exp(pi/2i)));
      if (min(min(min(dist0,dist1),dist2),dist3)==dist0)
         xbinstream_rx(1+2*(i-1))=0;
         xbinstream_rx(2+2*(i-1))=0;
      elseif (min(min(min(dist0,dist1),dist2),dist3)==dist1)
         xbinstream_rx(1+2*(i-1))=0;
         xbinstream_rx(2+2*(i-1))=1;
      elseif (min(min(min(dist0,dist1),dist2),dist3)==dist2)
         xbinstream_rx(1+2*(i-1))=1;
         xbinstream_rx(2+2*(i-1))=0;
      else
         xbinstream_rx(1+2*(i-1))=1;
         xbinstream_rx(2+2*(i-1))=1;
      end;
   end;

   BER_file=0;
   for i=1:1:length(x_bin_str)
      if (x_bin_str(i)~=xbinstream_rx(i))
         BER_file=BER_file+1;
      end;
   end;
   BER_file=BER_file/length(x_bin_str);

   for i=1:1:N
      ststruct_rx(i)=bi2de(xbinstream_rx(1+(i-1)*sym_bits:sym_bits+(i-1)*sym_bits),'left-msb');
   end;

   disp('---------------------------------------');
   disp(' Used SNR for file Tx ');
   disp(SNR);
   disp(' Achieved BER (%) ');
   disp(100*BER_file);
   
   SNR = 15;
   
      for i=1:1:length(x_sign_n_15)
      dist0=abs(x_sign_n_15(i)-amplitude*(1+exp(pi/2i)));
      dist1=abs(x_sign_n_15(i)-amplitude*(1-exp(pi/2i)));
      dist2=abs(x_sign_n_15(i)-amplitude*(-1+exp(pi/2i)));
      dist3=abs(x_sign_n_15(i)-amplitude*(-1-exp(pi/2i)));
      if (min(min(min(dist0,dist1),dist2),dist3)==dist0)
         xbinstream_rx_15(1+2*(i-1))=0;
         xbinstream_rx_15(2+2*(i-1))=0;
      elseif (min(min(min(dist0,dist1),dist2),dist3)==dist1)
         xbinstream_rx_15(1+2*(i-1))=0;
         xbinstream_rx_15(2+2*(i-1))=1;
      elseif (min(min(min(dist0,dist1),dist2),dist3)==dist2)
         xbinstream_rx_15(1+2*(i-1))=1;
         xbinstream_rx_15(2+2*(i-1))=0;
      else
         xbinstream_rx_15(1+2*(i-1))=1;
         xbinstream_rx_15(2+2*(i-1))=1;
      end;
   end;

   BER_file=0;
   for i=1:1:length(x_bin_str)
      if (x_bin_str(i)~=xbinstream_rx_15(i))
         BER_file=BER_file+1;
      end;
   end;
   BER_file=BER_file/length(x_bin_str);

   for i=1:1:N
      ststruct_rx_15(i)=bi2de(xbinstream_rx_15(1+(i-1)*sym_bits:sym_bits+(i-1)*sym_bits),'left-msb');
   end;

   disp('---------------------------------------');
   disp(' Used SNR for file Tx ');
   disp(SNR);
   disp(' Achieved BER (%) ');
   disp(100*BER_file);

%(vii) 

   out_string=native2unicode(ststruct_rx,'ASCII');

   filename_out='shannon_even_rx_5db.txt';
   fid=fopen(filename_out,'w');
   outf=fprintf(fid,'%c',out_string);
   ld=fclose(fid);

   disp('---------------------------------------');
   disp(' Recieved ASCII ');
   disp('---------------------------------------');
   disp(out_string)
   
   out_string=native2unicode(ststruct_rx_15,'ASCII');

   filename_out='shannon_even_rx_15db.txt';
   fid=fopen(filename_out,'w');
   outf=fprintf(fid,'%c',out_string);
   ld=fclose(fid);

   disp('---------------------------------------');
   disp(' Recieved ASCII ');
   disp('---------------------------------------');
   disp(out_string)