function spe(filename)



[y, fs]=audioread(filename);             % read wave file

left=y(:,1);                        % Left channel 
right=y(:,2); % Right channel

fy = fft(left);                     % transform waveform to frequency domain

figure;
spectrogram(fy)

end
