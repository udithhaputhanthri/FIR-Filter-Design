clc;
clear all;
close all;

A=2;B=0;C=8;

Ap=0.03 + (0.01 * A);
Aa=45+B;
Op1=(C * 100) + 400;
Op2=(C * 100) + 950;
Oa1=(C * 100) + 500;
Oa2=(C * 100) + 800;
Os=2*((C * 100) + 1300);
Ap
Aa
Bt=min([(Oa1-Op1),(Op2-Oa2)]);

Oc1=Op1+ Bt/2;
Oc2=Op2-Bt/2;
T=2*pi/Os;

%compute h[n]
%%%

delta_p= (10^(0.05*Ap)-1)/(10^(0.05*Ap)+1);
delta_a= 10^(-0.05*Aa);
delta=min(delta_p,delta_a);
Aa=-20*log10(delta);
Ap=20*log10((1+delta)/(1-delta));
if Aa<=21
    alpha=0;
    D=0.9222;
elseif 21<Aa<=50
    alpha=0.5842*(Aa-21)^0.4+0.07886*(Aa-21);
    D=(Aa-7.95)/14.36;
else
    alpha=0.1102*(Aa-8.7);
    D=(Aa-7.95)/14.36;
end

N=ceil(Os*D/Bt+1)+~mod(ceil(Os*D/Bt+1),2);

%%%%%%%%%%%%%%
n=[-(N-1)/2:1:(N-1)/2];
%%%%%%%%%%%%%%

beta=alpha*(1-((2.*n)./(N-1)).^2).^0.5;
h_d=(sin(Oc1*T*n)-sin(Oc2*T*n))./(n*pi);
h_d((N+1)/2)=1+2*(Oc1-Oc2)/Os;


inf_=100;
w= I_note(beta,inf_)./I_note(alpha,inf_);

h=w.*h_d;

figure;
stem(w);
title('Kaiser Window - Time Domain');
xlabel('n');
ylabel('Amplitude');
figure;
stem(h);
title('Kaiser Window - Filter Causal Impulse Response- Time Domain');
xlabel('n');
ylabel('Amplitude');
hold;
[H, H_freq_DT]=freqz(h,1);
H_=fft(h);
freq=H_freq_DT*Os/(2*pi);
figure;
H_mag_log=20*log10(abs(H));
plot(freq,H_mag_log);
xlabel('Frequency');
ylabel('Amplitude (dB)');
title('Kaiser Window - Filter Magnitude Response- Frequency domain');


%%%%%%%%%%%%%%%%%%%%%%%%%% 2----ANOTHER TWO PLOTS - LOW/HIGH ...
figure;
plot(freq,H_mag_log);
axis([800,1300,-0.1,0.1])
xlabel('Frequency');
ylabel('Amplitude (dB)');
title('Kaiser Window - LowerPassband- Magnitude Response');
figure;
plot(freq,H_mag_log);
axis([1600,2100,-0.1,0.1])
xlabel('Frequency');
ylabel('Amplitude (dB)');
title('Kaiser Window - UpperPassband- Magnitude Response');


w_rect=[];
%h_rect=w_rect.*h_d;
h_rect=h_d;
figure;
stem(h_rect);
title('Rectangular Window - Filter Causal Impulse Response- Time Domain');
xlabel('n');
ylabel('Amplitude');
hold;
[H_rect, H_freq_DT_rect]=freqz(h_rect,1);
freq_rect=H_freq_DT_rect*Os/(2*pi);
figure;
H_mag_log_rect=20*log10(abs(H_rect));
plot(freq_rect,H_mag_log_rect);
xlabel('Frequency');
ylabel('Amplitude (dB)');
title('Rectangular Window - Filter Magnitude Response- Frequency domain');



%%%%%%%%%%%%%%%%%%%%%%% INPUT - OUTPUT
O1=Oc1/2; 
O2=(Oc1+Oc2)/2;
O3=(Oc2+Os/2)/2;
width=[1:1000]; 
x=cos(O1*T*width)+cos(O2*T*width)+cos(O3*T*width);
x_ideal=cos(O1*T*width)+cos(O3*T*width);

figure;
[X,freq_DT]=freqz(x,1);
freq=freq_DT*Os/(2*pi);
plot(freq,abs(X));
xlabel('Frequency');
ylabel('Amplitude');
title('Input Signal- Frequency Domain');

figure;
stem(x);
xlabel('n');
ylabel('Amplitude');
title('Input Signal- Time Domain');

figure;
X_filtered=X.*H;
x_filtered=abs(ifft(X_filtered));
plot(freq,abs(X_filtered));
title('Output Signal- Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');

figure;
[X_ideal,freq_DT_ideal]=freqz(x_ideal,1);
freq_id=freq_DT_ideal*Os/(2*pi);
plot(freq_id,abs(X_ideal));
title('Ideal Output Signal- Frequency Domain');
xlabel('Frequency');
ylabel('Amplitude');



[u_,d_]=invfreqz(X_filtered,freq_DT,length(x)-1,1);
figure;
subplot(311);stem(x_ideal);
title('Ideal Output Signal- Time Domain');
xlabel('n');
ylabel('Amplitude');
subplot(312);stem(u_);
title('Output Signal- Time Domain');
xlabel('n');
ylabel('Amplitude');
axis([0 1000 -2 2]);




%%%% Impulse Response by fvtool
%%%% Freq represetation by fvtool
fvtool(h);

%%%% Impule Resp-designfilt()
%%%% Freq resp-designfilt()

bandstop = designfilt('bandstopfir', 'PassbandFrequency1', Op1, 'StopbandFrequency1', Oa1, 'StopbandFrequency2', Oa2, 'PassbandFrequency2', Op2, 'PassbandRipple1', Ap, 'StopbandAttenuation', Aa, 'PassbandRipple2', Ap, 'SampleRate', 4200, 'DesignMethod', 'kaiserwin');
fvtool(bandstop);





