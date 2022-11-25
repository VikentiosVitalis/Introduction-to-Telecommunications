% Vikentios Vitalis el18803
% A = 8 + 3 = 11 = 1 + 1 = 2

clc;
clear;

% Question a

   nbits=46;
   Tbit=0.5;

   randomseq=(1+sign(randn(nbits,1)))/2;

   amplitude=2;
   sampperbit=16;
   
   v=0;
   for i=1:1:nbits
      for j=1:1:sampperbit
         v=v+1;
         time_pt(v)=(Tbit/sampperbit)*(v-1);
         if (randomseq(i)==0)
            x_sign(v)=0;
         else
            x_sign(v)=amplitude;
         end;
      end;
   end;

 % Graphs

   figure(1)
   stem(randomseq,'o');
   grid;
   xlabel('Bit Number of Sequence');
   ylabel('Bit Value');
   title('Random Sequence');

   figure(2)
   plot(time_pt,x_sign,'LineWidth',1.2);
   grid;
   axis([0 nbits*Tbit+1 -1 amplitude+1]);
   xlabel('Time (s)');
   ylabel('Signal Value');
   title('B-PAM Modulated Random Sequence bits');

% Question b

   PAMst=2;	
   symamp(1)=0;
   symamp(2)=amplitude;
   
   for i=1:1:PAMst   
      xstar_nonoise(i)=(i-1)*(Tbit*symamp(i)^2/Tbit)^0.5;
      ystar_nonoise(i)=0;
   end;
   
   figure(3)
   plot(xstar_nonoise,ystar_nonoise,'o');
   grid;
   axis([-max(xstar_nonoise)-1 max(xstar_nonoise)+1 -max(xstar_nonoise)-1 max(xstar_nonoise)+1]);
   xlabel('(Eb/Tb)^0.5 (I)');
   ylabel('(Eb/Tb)^0.5 (Q)');
   title('Constellation Map for B-PAM');

% Question c

   SNR=5;

   x_sign_energy=sum(abs(x_sign).^2);
   
   noise_sign_energy=x_sign_energy/10^(SNR/10);
   noise_length=length(x_sign);
   sigma_n=(noise_sign_energy/noise_length)^0.5;
   noise_sigr=normrnd(0,sigma_n,[1,length(x_sign)]);
   noise_sign=noise_sigr;   
   noise_sign_energy_calc=sum(abs(noise_sign).^2);
   
   SNRcalculation=10*log10(x_sign_energy/noise_sign_energy_calc);
   SNRresult=SNR-SNRcalculation;
   noise_sign=10^(-SNRresult/10)*noise_sign;
   
   x_sign_n=x_sign+noise_sign;
   
   figure(4)
   plot(time_pt,x_sign_n,'LineWidth',1.2);
   grid;
   axis([0 nbits*Tbit+1 -1-max(x_sign_n) 1+max(x_sign_n)]);
   xlabel('Time (s)');
   ylabel('Signal Value');
   title('B-PAM Signal + AWGN (Eb/No=5 db)');

   SNR=15;

   x_sign_energy=sum(abs(x_sign).^2);
   
   noise_sign_energy=x_sign_energy/10^(SNR/10);
   noise_length=length(x_sign);
   sigma_n=(noise_sign_energy/noise_length)^0.5;
   noise_sigr=normrnd(0,sigma_n,[1,length(x_sign)]);
   noise_sign=noise_sigr;   
   noise_sign_energy_calc=sum(abs(noise_sign).^2);
   
   SNRcalculation=10*log10(x_sign_energy/noise_sign_energy_calc);
   SNRresult=SNR-SNRcalculation;
   noise_sign=10^(-SNRresult/10)*noise_sign;
   
   x_sign_n=x_sign+noise_sign;
   
   figure(5)
   plot(time_pt,x_sign_n,'LineWidth',1.2);
   grid;
   axis([0 nbits*Tbit+1 -1-max(x_sign_n) 1+max(x_sign_n)]);
   xlabel('Time (s)');
   ylabel('Signal Value');
   title('B-PAM Signal + AWGN (Eb/No=15 db)');

% Question d 

   SNR=5;

   x_sign_energy=sum(abs(x_sign).^2);
   
   noise_sign_energy=x_sign_energy/10^(SNR/10);
   noise_length=length(x_sign);
   sigma_n=(noise_sign_energy/noise_length)^0.5;
   noise_sigr=normrnd(0,sigma_n,[1,length(x_sign)]);
   noise_sigi=normrnd(0,sigma_n,[1,length(x_sign)]);
   noise_sign=noise_sigr+noise_sigi*exp(pi/2i);
   noise_sign_energy_calc=sum(abs(noise_sign).^2);
   
   SNRcalculation=10*log10(x_sign_energy/noise_sign_energy_calc);
   SNRresult=SNR-SNRcalculation;
   noise_sign=10^(-SNRresult/10)*noise_sign;
   
   x_sign_n=x_sign+noise_sign;

   for i=1:1:nbits
      bit_energyr(i)=sum(Tbit/sampperbit*real(x_sign_n(1+(i-1)*sampperbit:sampperbit+(i-1)*sampperbit)).^2);
      bit_energyi(i)=sum(Tbit/sampperbit*imag(x_sign_n(1+(i-1)*sampperbit:sampperbit+(i-1)*sampperbit)).^2);
   end;
   
   for i=1:1:nbits
      xstar(i)=(bit_energyr(i)/Tbit)^0.5;
      ystar(i)=(bit_energyi(i)/Tbit)^0.5;
   end;
   
   figure(6)
   plot(xstar,ystar,'x');
   hold on
   plot(xstar_nonoise,ystar_nonoise,'o');
   hold off
   grid;
   axis([-max(xstar)-1 max(xstar)+1 -max(xstar)-1 max(xstar)+1]);
   xlabel('(Eb/Tb)^0^.^5 (I)');
   ylabel('(Eb/Tb)^0^.^5 (Q)');
   title('Constellation Map for B-PAM (Eb/No=5 db)');
   legend('x=Stars for signal with noise','o=Stars for noiseless signal');
      
   SNR=15;

   x_sign_energy=sum(abs(x_sign).^2);
   
   noise_sign_energy=x_sign_energy/10^(SNR/10);
   noise_length=length(x_sign);
   sigma_n=(noise_sign_energy/noise_length)^0.5;
   noise_sigr=normrnd(0,sigma_n,[1,length(x_sign)]);
   noise_sigi=normrnd(0,sigma_n,[1,length(x_sign)]);
   noise_sign=noise_sigr+noise_sigi*exp(pi/2i);
   noise_sign_energy_calc=sum(abs(noise_sign).^2);
   
   SNRcalculation=10*log10(x_sign_energy/noise_sign_energy_calc);
   SNRresult=SNR-SNRcalculation;
   noise_sign=10^(-SNRresult/10)*noise_sign;
   
   x_sign_n=x_sign+noise_sign;

   for i=1:1:nbits
      bit_energyr(i)=sum(Tbit/sampperbit*real(x_sign_n(1+(i-1)*sampperbit:sampperbit+(i-1)*sampperbit)).^2);
      bit_energyi(i)=sum(Tbit/sampperbit*imag(x_sign_n(1+(i-1)*sampperbit:sampperbit+(i-1)*sampperbit)).^2);
   end;
   
   for i=1:1:nbits
      xstar(i)=(bit_energyr(i)/Tbit)^0.5;
      ystar(i)=(bit_energyi(i)/Tbit)^0.5;
   end;
   
   figure(7)
   plot(xstar,ystar,'x');
   hold on
   plot(xstar_nonoise,ystar_nonoise,'o');
   hold off
   grid;
   axis([-max(xstar)-1 max(xstar)+1 -max(xstar)-1 max(xstar)+1]);
   xlabel('(Eb/Tb)^0^.^5 (I)');
   ylabel('(Eb/Tb)^0^.^5 (Q)');
   title('Constellation Map for B-PAM (Eb/No=15 db)');
   legend('x=Stars for signal with noise','o=Stars for noiseless signal');
   
% Question e 
   
   format long;
   
   SNR_start=0;
   SNR_stop=15;
   SNR_step=1;
   SNR_iter=floor((SNR_stop-SNR_start)/SNR_step)+1;

   nbits=600000;
   Tbit=0.5;

   randomseq=(1+sign(randn(nbits,1)))/2;

   amplitude=9;
   sampperbit=1;
   
   v=0;
   for i=1:1:nbits
      for j=1:1:sampperbit
         v=v+1;
         time_pt(v)=(Tbit/sampperbit)*(v-1);
         if (randomseq(i)==0)
            x_sign(v)=0;
         else
            x_sign(v)=amplitude;
         end;
      end;
   end;
   
   x_sign_energy=sum(abs(x_sign).^2);
   
   for s=1:1:SNR_iter

      SNR_st(s)=SNR_start+(s-1)*SNR_step;
      SNR=SNR_st(s);
   
      noise_sign_energy=x_sign_energy/10^(SNR/10);
      noise_length=length(x_sign);
      sigma_n=(noise_sign_energy/noise_length)^0.5;
      noise_sigr=normrnd(0,sigma_n,[1,length(x_sign)]);
      noise_sigi=normrnd(0,sigma_n,[1,length(x_sign)]);
      noise_sign=noise_sigr+noise_sigi*exp(pi/2i);
      noise_sign_energy_calc=sum(abs(noise_sign).^2);
   
      SNRcalculation=10*log10(x_sign_energy/noise_sign_energy_calc);
      SNRresult=SNR-SNRcalculation;
      noise_sign=10^(-SNRresult/10)*noise_sign;
   
      x_sign_n=x_sign+noise_sign;

      for i=1:1:nbits
         if (abs(x_sign_n(i))<(amplitude/2))
            rand_seq_rx(i)=0;
         else
            rand_seq_rx(i)=1;
         end;
      end;
   
      BER(s)=0;
      for i=1:1:nbits
         if (randomseq(i)~=rand_seq_rx(i))
            BER(s)=BER(s)+1;
         end;
      end;
      BER_st(s)=BER(s)/nbits;
            
   end;
   
   figure(8)
   semilogy(SNR_st,BER_st,'LineWidth',2);
   grid;
   xlabel('SNR (db)');
   ylabel('BER');
   title('BER vs SNR for B-PAM');

