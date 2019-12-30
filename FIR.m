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

Bt=min([(Oa1-Op1),(Op2-Oa2)]);

Oc1=Op1+ Bt/2;
Oc2=Op2-Bt/2;
T=2*pi/Os;

%compute h[n]
%%%

delta_p= (10^(0.05*Ap)-1)/(10^(0.05*Ap)+1);
delta_a= 10^(-0.05*Ap);
delta=min([delta_p,delta_a]);

Aa=-20*log(delta);

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
n=[-(N-1)/2:(N-1)/2];
%%%%%%%%%%%%%%

beta=alpha*(1-((2.*n)./(N-1)).^2).^0.5;
h_d=(sin(Oc1*T*n)-sin(Oc2*T*n))./(n*pi);
h_d((N+1)/2)=1+2*(Oc1-Oc2)/Os;


inf_=100;
w= I_note(beta,inf_)./I_note(alpha,inf_);

h=w.*h_d;

m=1;m_=3;
figure;
subplot(m,m_,1);
stem(w);
title('window (Question_1)');
subplot(m,m_,2);
stem(h_d);
%plot(interp1(h_d,n));
title('h_d');
subplot(m,m_,3);
stem(h);
title('h (Question_2)');
hold;
figure;

%%%%%%%%%%%%%%%%%%%%%%%%%%%
freqz(h_d,1);
[H, H_freq_DT]=freqz(h,1);
freq=H_freq_DT*Os/(2*pi);
subplot(2,2,1);
plot(freq,abs(H));
xlabel('Frequency')
ylabel('Amplitude')
title('Kaiser Window - Freq Response');

subplot(2,2,2);
H_mag_log=20*log10(abs(H));
plot(freq,H_mag_log);
xlabel('Frequency')
ylabel('Amplitude (dB)')
title('Kaiser Window - Freq Response (Question_3,4)');



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
O1=Oc1/2;
O2=(Oc1+Oc2)/2;
O3=(Oc2+Os/2)/2;
k=[1:length(H)]; %%%%%%%%%%% CORRECT THISS !!!!!!!!!!!!
x=cos(O1*T*k)+cos(O2*T*k)+cos(O3*T*k);


%figure;
%subplot(1,2,1);
%stem(x);
figure;
subplot(1,2,1);
[X,freq_DT]=freqz(x,1);
freq=freq_DT*Os/(2*pi);
plot(freq,X)
xlabel('Frequency')
ylabel('Amplitude')
title('X');
subplot(1,2,2);
X_filtered=X.*H;
x_filtered=ifft(X_filtered);
plot(abs(X_filtered));
title('X filtered frequencyresponse');
xlabel('Frequency')
ylabel('Amplitude')
hold;





