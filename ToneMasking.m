%% Tone Masking With Noise
clear;
fs=16000;
samples = 10000;

%% Generate Gaussian Noise.

noise = wgn(samples,1,0); %white gaussian noise

%% Plot Gaussian Noise.

figure(1),plot(noise),title('White Gaussian Noise'), xlabel('Sample No.'), ylabel('amplitude')

%% Bandstop filter: Bandwidth Ajust 
% Adjust this variable to adjust the width of the Bandstop filter. Reduce
% the filter width until you can no longer hear the sinewave

BWadjust=2; %This variable is used to adjust the critical band.
% debug note. Min. acceptable value = 1.2

%% Bandstop filter construction

BWadjust=BWadjust/100; %Dividing by 100 is used because the inversely 
%proportional relationship of the variable BWadjust and the effects on Fp1
%and Fp2.

Fp1=0.5-BWadjust;     %Fp1=1/16;
Fp2=0.5+BWadjust;   %Fp2=15/16;
Fst1=Fp1+0.01;   % Fst1=2/16;
Fst2=Fp2-0.01;   % Fst2=14/16;

F1_LowCutoff=Fp1*(fs/2);
F2_HighCuttoff=Fp2*(fs/2);

Ap1=1;
Ast=60;
Ap2=1;

%Bandstop filter object using fdesign.
d = fdesign.bandstop('Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2', ... 
    Fp1,Fst1,Fst2,Fp2,Ap1,Ast,Ap2); 

%Final filter variable using function design().
Hd = design(d,'equiripple'); %FIR equiripple band stop filter 

FreqLow=(2*pi)/16;
FreqHi=(14*pi)/16;

%% Plot magnitude response of the bandstop filter. 
fvtool(Hd),

%% Calculate the stop band borders, center frequency and width in Hz 

freqHzLow=(fs*FreqLow)/(2*pi);
freqHzHigh=(fs*FreqHi)/(2*pi);
FilterWidth = freqHzHigh-freqHzLow;
CenterFreq = (freqHzLow+freqHzHigh)/2;

%% Apply the band stop filter to the Gaussian noise 

filtered_noise = filter(Hd,noise);

figure(3),subplot(1,2,1), plot(noise), xlabel('Sample No.'), ylabel('amplitude'), title('Un-filtered Gaussian Noise')
subplot(1,2,2), plot(filtered_noise), xlabel('Sample No.'), ylabel('amplitude'), title('Filtered Gaussian Noise')

%% Plot filtered noise.  

noisedft = fft(noise); 
filtered_noisedft = fft(filtered_noise); 
freq = 0:(2*pi)/length(noise):pi; 

figure(4),plot(freq/pi,abs(noisedft(1:length(noise)/2+1))),hold on, 
plot(freq/pi,abs(filtered_noisedft(1:length(filtered_noise)/2+1)),'r','linewidth',2),hold off, 
xlabel('Normalized Frequency (\times\pi rad/sample)'), ylabel('Magnitude'), legend('Original Signal','Bandstop Signal'),title('FFT spectra of the original and filtered noise signals')

%% Generate a sinusoid at the center of Bandpass filter frequencies. 

noise=noise'; %rearrange matrix to suit dimensions of y in preparation for addition.
amp=0.1*std(noise); %noise amplitude adjustment
ts=1/fs;
t=1:1:10000;
wCenter=(CenterFreq*2*pi)/fs;
y=sin(wCenter*t)*amp;

%% Play Sinusoid and Plot Sinusoid.
%sound(y,fs); %uncomment to listen.
figure(4),plot(y),xlim([0 1000]),xlabel('Sample No.'),ylabel('amplitude'),title('4kHz Sinewave')

%% Add Sine wave to gaussian noise.
output=y+noise;

%% Play and plot noisy sinusoid.

figure(5), plot(output), title('Sinewave with filtered Noise'), xlim([0 1000]), xlabel('Sample No.'), ylabel('amplitude')

sinusoiddft = fft(y);
maxsin=max(sinusoiddft);
filtered_noisedft = fft(filtered_noise);
freq2 = 0:(2*pi)/length(noise):pi;

figure(6),
plot(freq/pi,abs(sinusoiddft(1:length(noise)/2+1))), hold on
plot(freq/pi,abs(filtered_noisedft(1:length(filtered_noise)/2+1)),'r','linewidth',2),hold off, grid on,
xlabel('Normalized Frequency (\times\pi rad/sample)'), ylabel('Magnitude'), legend('Sinewave','Bandstop Signal'),title('Noise+4kHz Sinewave')

FiltNoiseSin=y+filtered_noise';

%Normalise audio output. Peak value 0 dBFS.
pkch=max(FiltNoiseSin);
norm=FiltNoiseSin*(1/pkch);

%% Sound output and info
sound(norm,fs); %uncomment to play audio. Keep speaker volume low before playing. Take note of BP filter width when sin note can no longer be heard.
filterwidth=F2_HighCuttoff-F1_LowCutoff;

fprintf('Sine tone freq: %d Hz\n', CenterFreq);
fprintf('BP Filter Width: %d Hz\n', filterwidth);
