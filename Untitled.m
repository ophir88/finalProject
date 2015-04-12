



fs = 2000;
[y, fs]=audioread('media/music notes/center_do_re_mi.wav'); % Read wave file
left=y(:,1); % Left channel 
right=y(:,2); % Right channel
spectrogram(left, ones(1, 256),0,256, 2000, 'yaxis')

[funky,fs] = wavread('media/music notes/center_do_re_mi.wav');
wavplay(funky, fs);