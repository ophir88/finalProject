function main2(filename, midiname)


% ---------------- variables: ---------------
midi = readmidi(midiname);
Notes = midiInfo(midi,0);
tempo = midi.ticks_per_quarter_note;
bpm =  60*1000/tempo;
bpm = 96;
% read file:
[mySignal,Fs] = audioread(filename);
% sound(mySignal,Fs)
% mySignal = mySignal/max(abs(mySignal));

% file length in seconds:
signalTime = length(mySignal)/Fs;
% file length vector
L = length(mySignal);

%window size in seconds: (set according to midi file)

windowSizeTime = 60/bpm; %seconds
% windowSizeTime = 0.5 %seconds
%window size vector:
windowSize = round(windowSizeTime*Fs);

%   overlap: currently not used:
% overlapTime = 0.01;
% overlapSize = overlapTime*Fs;

% threshold
threshold = .17;



% results: [between time, to time, note, amp]
results = [0,0,0,0,0,0,0,0];
results1 = zeros(floor(L/windowSize)+1,8);
ok = zeros(floor(L/windowSize)+1,3);
okIdx =1;
% keep temp peaks during 'for' loop: [between time, and time,freq, amp]
peaks = zeros(6,8);

peaksIdx = 1;
resultIdx = 1;
newFreq = 0;
notesIdx = 1;
notesSegIdx = 0;
desiredNotes = zeros(6,1);
prevDesiredNotes = zeros(6,1);

historyEnv = 0;
% ------------- start of function: ---------------------

% break the input into 'window size' elements, and analyze separatly
%   iterate over input in 'window size' jumps
for k=1 : windowSize : L
    desiredNotes = zeros(6,1);
    theEnd=0;
    ContinousNote =0;
    %     if (k+windowSize>L)
    %         partSignal = mySignal(k:L)
    %         theEnd = 1;
    %     else
    %         partSignal = mySignal(k:k+windowSize);
    %     end
    
    
    %set note we look for:
    j=1;
    for i=notesIdx:length(Notes)
        
        if j<7
            if (round(Notes(notesIdx,5)/windowSizeTime)==notesSegIdx)
                desiredNotes(j)=Notes(notesIdx,3);
                notesIdx=notesIdx+1;
                j=j+1;
                
            else
                j=7;
                break;
            end
        else
            break
        end
        
        
    end
    notesSegIdx=notesSegIdx+1;
    desiredNotes;
    
    
    % % % % %     % --*-*-**-**-*-*-*--* find envelope *-*-*-*-*-*-*-*

    % Normalize y; that is, scale all values to its maximum. Note how simple it is
    % to do this in MatLab.
    sound(mySignal(k:k+windowSize),Fs)
    y = abs(mySignal(k:k+windowSize));
    %         y = abs(partSignal);
  
   
    maxenv = y;
    
  
    peakwindowsize = 512;
    
    
    % Go through the input signal, taking the maximum value of the previous
    % peakwindowsize number of samples. We do this in the same way as we did the
    % average, but now we use max instead of sum.
    
    for p = peakwindowsize:length(maxenv)
        maxenv(p) = max(y(p-peakwindowsize+1:p));
    end
    
    
  
    
    % find amp peaks:
    peaksLoc = zeros(4,2);
    inPeak = 0;
    peakStartLoc = 0;
%     peakEnd = 0;
    peakEndLoc= 0;
    thresh = 0.38;
    peakInd =1;
    continuousNote = 0;
    for i=1:round(windowSize)       
        if (inPeak == 0)
            if (maxenv(i)>thresh)
                inPeak = 1;
                peakStartLoc = i;
            end
        else
            if (maxenv(i)<thresh)
%                 peakEnd = 1;
                inPeak = 0;
                peakEndLoc = i;
                if(peakEndLoc-peakStartLoc>3000)
                    
                    peaksLoc(peakInd,1)=peakStartLoc;
                    peaksLoc(peakInd,2)=peakEndLoc;
                    peakInd = peakInd+1;
                end
            end
        end 
    end
    if(inPeak == 1 && peakStartLoc<25000)
        peaksLoc(peakInd,1)=peakStartLoc;
        peaksLoc(peakInd,2)=windowSize;
        continuousNote = 1;
    end
    


      figure(round(k/windowSize)+101)
    plot (maxenv, 'm')
    
    for i=1:peakInd
    if (peaksLoc(i,2)-peaksLoc(i,1)>0)
        hold on
        plot([peaksLoc(i,1) peaksLoc(i,2)],[thresh thresh],'--or')

%         hline = gline; % Connect circles
%         set(hline,'Color','r')
    else 
        break
    end
    end
    
        's';

    

    
    
    
    
    % loop over peak locations, looking for frequencies:
    
    for i=1:peakInd
        
    if (peaksLoc(i,2)-peaksLoc(i,1)>0)
        
        
        
    else 
        break
    end
    end    
    
    
    
    
    NFFT = 2^nextpow2(windowSize); % Next power of 2 from length of y
    % get segment (check if not last one)
    if (k+windowSize>L)
        spectrum = fft(mySignal(k:L),NFFT);
    else
        spectrum = fft(mySignal(k:k+windowSize),NFFT);
    end
    spectrum = abs(spectrum(1:ceil(length(spectrum)/2)));
    
    % max(spectrum);
    filter = gaussfilt(128);
    spectrum = spectrum .^ 2;
    %     spectrum2 = conv(spectrum,filter);
    % regular fft:
    spectrum = spectrum/max(spectrum);
    % filtered fft:
    %     spectrum2 = spectrum2/max(spectrum2);
    % convert between vector location to real frequency:
    tranCoef = Fs/(length(spectrum)*2);
    %plot signals:
    %         figure(round(k/windowSize)+1)
    %                 figure(1)
    %
    %         plot(linspace(0,Fs/2,length(spectrum)),spectrum)
    %         axis([0 1500 0 1])
    %     figure(round(k/windowSize)+11)
    %     plot(linspace(0,Fs/2,length(spectrum2)),spectrum2)
    %     axis([0 1500 0 1])
    
    % ------------------- find peaks: ----------------
    
    % look only at relevant frequencies (piano notes)
    'loop';
    for i=round(27*length(spectrum)/Fs):round(4180*length(spectrum)/Fs)
        
        % max point:
        if(spectrum(i)>threshold & spectrum(i)>spectrum(i-1)& spectrum(i)>spectrum(i+1))
            
            if (peaksIdx>6)
                break;
            end
            % make sure peaks are not too close
            
            if(abs(i*tranCoef -peaks(peaksIdx,4)) > 10)
                
                % close peaks, look for stronger:
                % case new one is stronger, overwrite
                % older
                %                 if(spectrum(i) > peaks(peaksIdx,4))
                %new peak is close and much stronger, reaplace it with the
                %former one:
                newFreq = i*tranCoef;
                
                %check that it is not an harmonic of a already added note:
                octave = 0;
                for j=1:peaksIdx
                    
                    % if the new frequency is a multiplication of an
                    % existing one
                    
                    if(peaks(j,4)>0 &&(mod(newFreq,peaks(j,4)) < .03*newFreq || ...
                            peaks(j,4) - mod(newFreq,peaks(j,4))  < .03*newFreq))
                        
                        
                        %found octave, for now ignore:
                        octave = 1;
                        
                        % If the magnitude of the peak is strong enough, declare that it is a different note instead of a loud harmonic
                        if(spectrum(i) > .3*peaks(j,5))
                            
                            % add peak:
                            %                                                                 'before add'
                            %                                                                 peaks
                            %                                                                 'after add'
                            %start time
                            peaks(peaksIdx,1)=round(k/windowSize)*windowSizeTime;
                            %end time
                            peaks(peaksIdx,2)=(round(k/windowSize)+1)*windowSizeTime;
                            %segment number:
                            peaks(peaksIdx,3)=(round(k/windowSize));
                            %freq
                            peaks(peaksIdx,4)=i*tranCoef;
                            %ampli
                            peaks(peaksIdx,5)=spectrum(i);
                            % is harmony?
                            peaks(peaksIdx,6)=1;
                            % average amp
%                             peaks(peaksIdx,8)=average2;
                            
                            peaksIdx = peaksIdx +1;
                            
                            %                             'added'
                            break;
                        end
                        
                    else
                    end
                    
                    
                end
                if( octave == 0)
                    
                    'no octave';
                    peaks(peaksIdx,1)=round(k/windowSize)*windowSizeTime;
                    %end time
                    peaks(peaksIdx,2)=(round(k/windowSize)+1)*windowSizeTime;
                    %segment number:
                    peaks(peaksIdx,3)=(round(k/windowSize));
                    %freq
                    peaks(peaksIdx,4)=i*tranCoef;
                    %ampli
                    peaks(peaksIdx,5)=spectrum(i);
                    % is harmony?
                    peaks(peaksIdx,6)=0;
                    % average amp
%                     peaks(peaksIdx,8)=average2;
                    
                    
                    peaksIdx = peaksIdx +1;
                    if (peaksIdx==7)
                        break
                    end
                end
            else continue;
                %                 end
            end
            
        end
    end
    % ------------------- end find peaks -------------
    
    
    %     (horzcat(round(k/windowSize)*windowSizeTime,round((k+1)/windowSize)*windowSizeTime,peaks))
    
    
    
    
    % look for the desired notes in the notes we found:
    for l=1:length(desiredNotes)
        found = 0;
        if(desiredNotes(l)>0)
            for m=1:size(peaks)
                
                % if we found the desired note:
                if (desiredNotes(l)+3>peaks(m,4) &&desiredNotes(l)-3<peaks(m,4))
                    %segment number
                    ok(okIdx,1)=notesSegIdx;
                    %frequency
                    ok(okIdx,2)=desiredNotes(l);
                    % boolean 'was found':
                    ok(okIdx,3)=1;
                    okIdx=okIdx+1;
                    found = 1;
                    % note in the peaks matrix that this note is used:
                    peaks(m,7)=1;
                    break;
                end
                
            end
            % if note wasn't found:
            if (found==0)
                ok(okIdx,1)=notesSegIdx;
                ok(okIdx,2)=desiredNotes(l);
                ok(okIdx,3)=0;
                okIdx=okIdx+1;
            end
        else
            break
        end
    end
    
    % insert our found notes into a database: (look for unnecessary notes played):
    for m=1:size(peaks)
        if (peaks(m,5)>0)
            %start time
            results1(resultIdx,1)=peaks(m,1);
            %end time
            results1(resultIdx,2)=peaks(m,2);
            %segment number:
            results1(resultIdx,3)=peaks(m,3);
            %freq
            results1(resultIdx,4)=peaks(m,4);
            %amplitude
            results1(resultIdx,5)=peaks(m,5);
            % is a second harmony?
            results1(resultIdx,6)=peaks(m,6);
            % was it in the 'desired notes'
            % array?
            results1(resultIdx,7)=peaks(m,7);
            %average amplitude:
%             results1(resultIdx,8)=peaks(m,8);
            
            % increment index:
            resultIdx=resultIdx+1;
        else
            break
        end
        
    end
    
    peaks = zeros(6,8);
    peaksIdx = 1;
    prevDesiredNotes = desiredNotes;
    
    
    
end

[rows,columns] = size(results1);
if(columns >= 2)
    results1 = results1(:,1:columns);
end
[rows,columns] = size(ok);
if(columns >= 2)
    ok = ok(:,1:columns);
end
results1
ok
'end'

end


function out=fconv(f,h)
P=length(f)+length(h)-1;
fnew=zeros(1,P);
hnew=zeros(1,P);
fnew(1:length(f))=f;
hnew(1:length(h))=h;
out=ifft(fft(fnew).*fft(hnew));
end


function filter=gaussfilt(N)

% calculate alpha so the filter fills all N points
alpha=N;
first=-(1-N/2)*exp(-(1-N/2)^2/alpha);
count=0;
while first<.1*(-(1530/4000*N-N/2)*exp(-(1530/4000*N-N/2)^2/alpha))
    count=count+1;
    alpha=N*500*count;
    first=-(1-N/2)*exp(-(1-N/2)^2/alpha);
end


for n=1:N
    filter(n)=-(n-N/2)*exp(-(n-N/2)^2/alpha);   % d/dt of a gaussian
end
filter=filter/sum(abs(filter));     % normalization
end




