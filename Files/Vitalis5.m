% Vikentios Vitalis el18803
% 8 + 3 = 11 = 1 + 1 = 2

clc;
clear;
format short;

% 	Question a - Wav Read


   AM_sum=8+0+3;
   
   if (mod(AM_sum,2)==0)
      filename='soundfile2_lab2.wav';
   else
      filename='soundfile1_lab2.wav';
   end;   

   [st_struct,fs]=audioread(filename);   

   Amax=max(st_struct);
   Amin=min(st_struct);
   
   for i=1:1:length(st_struct)
      st_struct_norm(i)=floor(255*(st_struct(i)-Amin)/(Amax-Amin));
   end;

   sym_bits=8;
   for i=1:1:length(st_struct_norm)
      binfile(i,:)=de2bi(st_struct_norm(i),sym_bits,'left-msb');
      if (mod(i,1000)==0)
         clc;
         disp(' Binary Conversion Progress (%)');
         disp(100*i/length(st_struct_norm));
      end;
   end;
   clc;
   
   v=0;
   for i=1:1:length(st_struct)
      for j=1:1:sym_bits
         v=v+1;
         bin_stream(v)=binfile(i,j);
      end;
   end;

   figure(1)
   i=1:1:length(st_struct);
   semilogy((i-1)/fs,st_struct,'LineWidth',1);
   grid;
   xlabel(' Time (sec) ');
   ylabel(' Signal Value ');
   title(' Wav Signal Values');
   
   sound(st_struct,fs);


   % Question b - Quantizer


   Amax_norm=max(st_struct_norm);
   Amin_norm=min(st_struct_norm);
   n=8;
   N=length(st_struct_norm);

   delta=(Amax_norm-Amin_norm)/2^n;

   partition(1) = Amin_norm+delta/2;
   for i=1:1:2^n-1
      partition(i+1)=partition(i)+delta;
   end;

   x_ind=quantiz(st_struct_norm,partition);

   for i=1:1:N
      if (x_ind(i)+1>2^n)
         x_qnd(i)=partition(2^n)+delta/2;
      else
         x_qnd(i)=partition(x_ind(i)+1);
      end;
   end;

   sym_bits=8;
   x_bin=de2bi(floor(x_qnd),sym_bits,'left-msb');

   v=0;
   for i=1:1:length(st_struct)
      for j=1:1:sym_bits
         v=v+1;
         x_bin_stream(v)=x_bin(i,j);
      end;
   end;

   figure(2)
   i=1:1:length(x_qnd);
   semilogy((i-1)/fs,x_qnd,'LineWidth',1);
   grid;
   xlabel(' Time (sec) ');
   ylabel(' Signal Value ');
   title('Quantized of Normalized Wav Signal Values');


% 	Questions c-d - Q-PSK


   SNR=4;
   Amp=1;

   L=4;
   % Mapping vector for M - PSK Gray encoding
   k=log2(L); 			% Number of bits per point

   ph1=[pi/4]; 
   theta=[ph1; -ph1; pi-ph1; -pi+ph1]; 
   mapping=exp(1j*theta); 
   if(k>2) 
      for j=3:k 
         theta=theta/2; 
         mapping=exp(1j*theta); 
         mapping=[mapping; -conj(mapping)]; 
         theta=real(log(mapping)/1j); 
      end;
   end;
   xstar_noiseless=(Amp^2/2)^0.5*real(mapping);
   ystar_noiseless=(Amp^2/2)^0.5*imag(mapping);

   xsym=bi2de(reshape(x_bin_stream,log2(L),length(x_bin_stream)/log2(L)).','left-msb'); 
   x_map=[]; 
   for k=1:length(xsym) 
      x_val(k)=real(mapping(xsym(k)+1));
      y_val(k)=imag(mapping(xsym(k)+1));
      x_map(k)=x_val(k)+exp(pi/2i)*y_val(k); 
   end;
   x_map=(Amp^2/2)^0.5*x_map;

   x_sig_energy=sum(abs(x_map).^2);
   
   noise_sig_energy=x_sig_energy/10^(SNR/10);
   noise_length=length(x_map);
   sigma_n=(noise_sig_energy/noise_length)^0.5;
   noise_sigr=normrnd(0,sigma_n,[1,length(x_map)]);
   noise_sigi=normrnd(0,sigma_n,[1,length(x_map)]);
   noise_sig=noise_sigr+noise_sigi*exp(pi/2i);
   noise_sig_energy_calc=sum(abs(noise_sig).^2);
   
   SNR_calc=10*log10(x_sig_energy/noise_sig_energy_calc);
   SNR_res=SNR-SNR_calc;
   noise_sig=10^(-SNR_res/10)*noise_sig;
   
   x_sig_n=x_map+noise_sig;

   xstar=real(x_sig_n);
   ystar=imag(x_sig_n);
   xstar_noiseless=real(mapping);
   ystar_noiseless=imag(mapping);
   
   figure(3)
   plot(xstar,ystar,'x');
   hold on
   plot(xstar_noiseless,ystar_noiseless,'o');
   hold off
   grid;
   axis([-max(xstar)-1 max(xstar)+1 -max(xstar)-1 max(xstar)+1]);
   xlabel(' (Eb/Tb)^0.5 (I) ');
   ylabel(' (Eb/Tb)^0.5 (Q) ');
   title('Constellation Map for Q-PSK (Eb/No=14 db)');
   legend('x=Stars for signal with noise','o=Stars for noiseless signal');


% 	Question e - Demodulation of Q-PSK


   for i=1:1:length(x_sig_n)
      dist0=abs(x_sig_n(i)-Amp*(1+exp(pi/2i)));
      dist1=abs(x_sig_n(i)-Amp*(1-exp(pi/2i)));
      dist2=abs(x_sig_n(i)-Amp*(-1+exp(pi/2i)));
      dist3=abs(x_sig_n(i)-Amp*(-1-exp(pi/2i)));
      if (min(min(min(dist0,dist1),dist2),dist3)==dist0)
         x_bin_stream_rx(1+2*(i-1))=0;
         x_bin_stream_rx(2+2*(i-1))=0;
      elseif (min(min(min(dist0,dist1),dist2),dist3)==dist1)
         x_bin_stream_rx(1+2*(i-1))=0;
         x_bin_stream_rx(2+2*(i-1))=1;
      elseif (min(min(min(dist0,dist1),dist2),dist3)==dist2)
         x_bin_stream_rx(1+2*(i-1))=1;
         x_bin_stream_rx(2+2*(i-1))=0;
      else
         x_bin_stream_rx(1+2*(i-1))=1;
         x_bin_stream_rx(2+2*(i-1))=1;
      end;
   end;

% 	Question st - BER Calculation


   BER_file=0;
   for i=1:1:length(x_bin_stream)
      if (x_bin_stream(i)~=x_bin_stream_rx(i))
         BER_file=BER_file+1;
      end;
   end;
   BER_file=BER_file/length(x_bin_stream);

   for i=1:1:N
      st_struct_rx(i)=bi2de(x_bin_stream_rx(1+(i-1)*sym_bits:sym_bits+(i-1)*sym_bits),'left-msb');
   end;

   disp('---------------------------------------');
   disp(' Used SNR for file Tx ');
   disp(SNR);
   disp(' Achieved BER (%) ');
   disp(100*BER_file);
   
      

%  Question z - Write Rxc Wav to file

   out_values=(Amax-Amin)*st_struct_rx/255;

   if (mod(AM_sum,2)==0)
      filename_out='soundfile2_lab2_rx_14db.wav';
   else
      filename_out='soundfile1_lab2_rx_14db.wav';
   end;   

   figure(4)
   i=1:1:length(out_values);
   semilogy((i-1)/fs,out_values,'LineWidth',1);
   grid;
   xlabel(' Time (sec) ');
   ylabel(' Signal Value ');
   title(' Rx Wav Signal Values');

   disp('---------------------------------------');
   disp(' Original Wav file Info ');
   disp('---------------------------------------');
   audioinfo(filename)

   audiowrite(filename_out,out_values,fs,'BitsPerSample',8);
   disp('---------------------------------------');
   disp(' Rx Wav file Store Info ');
   disp('---------------------------------------');
   audioinfo(filename_out)
   
   sound(out_values,fs);   
   
