function main(filename, midiname, guitar, plotOn)


close all;
% ---------------- variables: ---------------
midi = readmidi(midiname);
Notes = midiInfo(midi,0);
NoteCount = size(Notes);
NoteCount = NoteCount(1);

bpm = midi.bpm;
bmp2 = 60000000/midi.micPq;
eigth = bmp2*2;



% read file:
[mySignal,Fs] = audioread(filename);


% file length in seconds:
signalTime = length(mySignal)/Fs;

% file length (vector)
L = length(mySignal);

%window size in seconds: (length of quarter note (according to midi file)
windowSizeTime = 60/bpm; %seconds

%window size vector:
windowSize = round(windowSizeTime*Fs);

% frequency amplitude threshold
threshold = .01;


% results: [time start, time end , segment number, frequency, amplitude, harmony]
results1 = zeros(floor(L/windowSize)+1,8);
% note comparison: what was played VS what needed to be played
ok = zeros(floor(L/windowSize)+1,3);
okIdx =1;

% temporary frequency array for calculations (is merged to 'results' at end of each loop):
peaks = zeros(6,8);

peaksIdx = 1;
resultIdx = 1;
newFreq = 0;
notesIdx = 1;
notesSegIdx = 0;
% listed of desired notes for each segment according to midi file:
desiredNotes = zeros(6,2);
% list of previous desired notes (will server for notes which last longer
% than 1 segemnt, i.e. half notes and up):
prevDesiredNotes = zeros(6,2);
historyEnv = 0;

%  variables to destenguish place in time:
noteIdx = 1;
first = 1;
NxtNoteTime = Notes(noteIdx,5);
currTime=0;
yEnd = 0;
pastTime=0;


% figure number:
figureNum = 1;
% ------------- Main function Start: ---------------------

% break the input into 'window size' elements, and analyze separatly:
% iterate over input in 'window size' jumps

% for k=1 : windowSize : L
while currTime*Fs<L
    %     currTime*Fs
    %     L
    
    %     ============== get desired notes for this iteration: ================
    
    tempDesiredNotes = zeros(6,1);
    tempDesiredIdx = 1;
    
    % look for all notes which start at the current time:
    if(noteIdx>=NoteCount)
        break
    end
    
    while Notes(noteIdx,5)<=currTime && tempDesiredIdx<=6
        % frequency:
        tempDesiredNotes(tempDesiredIdx,1) = Notes(noteIdx,3);
        % note off time:
        tempDesiredNotes(tempDesiredIdx,2) = Notes(noteIdx,6);
        noteIdx = noteIdx+1;
        tempDesiredIdx = tempDesiredIdx+1;
        if(tempDesiredIdx>6)
            break
        end
    end
    
    % look for all notes from previous desiredNotes, which are still
    % played in this window
    for l = 1 : size(desiredNotes)
        
        % if note is played at least 1/16 note into current window, add it
        if(desiredNotes(l,1)>0 && desiredNotes(l,2)-currTime>windowSize/(Fs*4) && tempDesiredIdx<=6)
            tempDesiredNotes(tempDesiredIdx,1) = desiredNotes(l,1);
            tempDesiredNotes(tempDesiredIdx,2) = desiredNotes(l,2);
            tempDesiredIdx = tempDesiredIdx+1;
        end
        if(tempDesiredIdx>6)
            break
        end
    end
    
    % look for all notes start close to the curr time (within 1/16 note away):
    %     a = windowSize/(Fs*2);
    while Notes(noteIdx,5)-currTime<= windowSize/(Fs*4) && tempDesiredIdx<=6
        % frequency:
        tempDesiredNotes(tempDesiredIdx,1) = Notes(noteIdx,3);
        % note off time:
        tempDesiredNotes(tempDesiredIdx,2) = Notes(noteIdx,6);
        noteIdx = noteIdx+1;
        tempDesiredIdx = tempDesiredIdx+1;
        
    end
    
    %     added all relevant notes, now update the window size we will be
    %     working on :
    %  !*!*!!*!*!*!* TODO (make sure to update currTime at the end!)
    
    
    %   if a new note is starting before the end of this window, set the end of
    %   the window to be the begining of the next note:
    if((currTime+windowSize/Fs) > Notes(noteIdx,5)+0.001)
        yEnd = Notes(noteIdx,5)*Fs;
    else
        %         o.w. set the end time to be current time + window size:
        yEnd = currTime*Fs+windowSize;
    end
    if (yEnd > L)
        yEnd = L;
    end
    desiredNotes = tempDesiredNotes;
    %     ================== setting variable for current loop ================
    
    %   setting currTime to be greater than zero (for accessing the signal)
    if(currTime==0)
        currTime = 0.0001;
    end
    yStart = round(currTime*Fs);
    %     sound the current window:
    sound(mySignal(yStart:yEnd),Fs)
    % before transformation to absolute value, create a copy for the
    % fourier analysis:
    yFourier = mySignal(yStart:yEnd);
    % get absolute value
    y = abs(yFourier);
    
    %     ====================== get envelop: ====================
    %   copy of signal window:
    maxenv = y;
    %   windowsize to look for maximum value at:
    peakwindowsize = 512;
    % get envelop vector:
    average = 0;
    's';
    % get envelop:
    for p = 1:length(maxenv)-peakwindowsize
        average = average+y(p);
        maxenv(p) = max(y(p:p+peakwindowsize));
    end
    if(length(maxenv)-peakwindowsize<=peakwindowsize)
        peakwindowsize=length(maxenv)-peakwindowsize-1;
    end
    if(length(maxenv)-peakwindowsize>1)
        for p = length(maxenv)-peakwindowsize:length(maxenv)
            %         p-peakwindowsize+1
            try
                maxenv(p) = max(y(p-peakwindowsize+1:p));
            catch exception
                's';
            end
            
            average = average+y(p);
        end
    end
    
    %     average value:
    average = average/length(y);
    
    %   ====================== find amplitue peaks: ====================
    %    ----- variables: ------
    %   set memory for 4 peaks:
    peaksLoc = zeros(4,2);
    inPeak = 0;
    peakStartLoc = 0;
    thresh = average;
    peakInd =1;
    %   if at the end of the window the peak is still high, set the note as
    %   continuos:
    continuousNote = 0;
    for i=1:length(maxenv)
        %         if not in peak and higher than thresh, set to "in peak" and save
        %         location of begining of peak
        if (inPeak == 0)
            if (maxenv(i)>thresh)
                inPeak = 1;
                peakStartLoc = i;
            end
        else
            %             if already  in peak, and reached the end of peak, set inpeak
            %             to false and save location of end of peak
            if (maxenv(i)<thresh)
                inPeak = 0;
                peakEndLoc = i;
                %                 if peak is at least the size of an eight note, save it
                if(peakEndLoc-peakStartLoc>windowSize/8)
                    peaksLoc(peakInd,1)=peakStartLoc;
                    peaksLoc(peakInd,2)=peakEndLoc;
                    peakInd = peakInd+1;
                end
            end
        end
    end
    %     if we finish the window and are still in peak, set as continuos note:
    if(inPeak == 1 && peakStartLoc<length(y)*7/8)
        peaksLoc(peakInd,1)=peakStartLoc;
        peaksLoc(peakInd,2)=length(y);
        continuousNote = 1;
    end
    
    
    %   ====================== envelop and peak ploting: ====================
    
    if(plotOn == 1)
        %     figure(figureNum+101);
        figure;
        subplot(1,2,1);
        plot (maxenv, 'm')
        for i=1:peakInd
            if (peaksLoc(i,2)-peaksLoc(i,1)>0)
                hold on
                plot([peaksLoc(i,1) peaksLoc(i,2)],[thresh thresh],'--or')
            else
                break
            end
        end
        hold off;
    end
    %   ====================== Frequency domain analysis: ====================
    
    % loop over peak locations, looking for frequencies:
    for j=1:peakInd
        %         make sure ther is a peak (i.e. it is not the default zeros)
        if (peaksLoc(j,2)-peaksLoc(j,1)>0)
            %  fourier transform:
            
            NFFT = 2^nextpow2(peaksLoc(j,2)-peaksLoc(j,1)); % Next power of 2 from length of y
            spectrum = fft(yFourier(peaksLoc(j,1):peaksLoc(j,2)),NFFT);
            spectrum = abs(spectrum(1:ceil(length(spectrum)/2)));
            spectrum = spectrum .^ 2;
            %           normalize:
            spectrum = spectrum/max(spectrum);
%             for i=1:length(spectrum)
%                 if(spectrum(i)<0.0001)
%                     spectrum(i)=0;
%                 end
%                 
%             end
            
%             downSpectrum1 = [spectrum(1:2:end) zeros(1,length(spectrum)/2)];
%             downSpectrum2 = spectrum(1:3:end);
%             downSpectrum2 = [downSpectrum2 zeros(1,length(spectrum)-length(downSpectrum2))];
%             finalSpec = spectrum.*downSpectrum1.*downSpectrum2;
%             figure;
%             subplot(2,2,1);
%             axis([0 1500 0 1]);
% 
%             plot(linspace(0,Fs/2,length(spectrum)),spectrum);
%             subplot(2,2,2);
%             plot(linspace(0,Fs/2,length(spectrum)),downSpectrum1);
%             subplot(2,2,3);
%             plot(linspace(0,Fs/2,length(spectrum)),downSpectrum2);
%             subplot(2,2,4);
%             plot(linspace(0,Fs/2,length(finalSpec)),finalSpec);
            
            %             tansform coefficient:
            tranCoef = Fs/(length(spectrum)*2);
            
            
            % ------------------- find peaks: ----------------
            % look only at relevant frequencies (note frequency spectrum)
            for i=round(27*length(spectrum)/Fs):round(4180*length(spectrum)/Fs)
                %variable for close peaks:
                minorPeak = 0;
                % max point:
                if(spectrum(i)>threshold && spectrum(i)>spectrum(i-1)&& spectrum(i)>spectrum(i+1))
                    
                    if (peaksIdx>6)
                        break;
                    end
                    % ------------------- check for close range peaks: ----------------
                    %if this is not the first peak
                    if(peaksIdx>1)
                        % check that this is the same peak window:
                        if(peaks(peaksIdx-1,8)==j)
                            %   if the location of this peak is close to the
                            %  location of the previous peak
                            if (abs((i-1)*tranCoef -peaks(peaksIdx-1,4)) < 15)
                                % if the amplitude of this peak is higher then
                                % previous peak, overwrite previous:
                                if (spectrum(i)>peaks(peaksIdx-1,5))
                                    peaksIdx = peaksIdx-1;
                                else
                                    minorPeak = 1;
                                end
                            end
                        end
                        
                        
                    end
                    
                    if(minorPeak == 0)
                        % close peaks, look for stronger:
                        % case new one is stronger, overwrite
                        % older
                        % if(spectrum(i) > peaks(peaksIdx,4))
                        %new peak is close and much stronger, reaplace it with the
                        %former one:
                        newFreq = (i-1)*tranCoef;
                        
                        % ------------------- check for harmonics: ----------------
                        
                        %check that it is not an harmonic of a already added note:
                        octave = 0;
                        for k=1:peaksIdx
                            % if the new frequency is a multiplication of an
                            % existing one
                            if(peaks(k,4)>0 &&(mod(newFreq,peaks(k,4)) < .03*newFreq || ...
                                    peaks(k,4) - mod(newFreq,peaks(k,4))  < .03*newFreq))
                                %found octave, for now ignore:
                                octave = 1;
                                % If the magnitude of the peak is strong enough, declare that it is a different note instead of a loud harmonic
                                if(spectrum(i) > .3*peaks(k,5))
                                    %start time
                                    peaks(peaksIdx,1)=currTime +  peaksLoc(j,1)/Fs;
                                    %end time
                                    peaks(peaksIdx,2)=currTime + peaksLoc(j,2)/Fs;
                                    %segment number:
                                    peaks(peaksIdx,3)=figureNum;
                                    %freq
                                    peaks(peaksIdx,4)=newFreq;
                                    %ampli
                                    peaks(peaksIdx,5)=spectrum(i);
                                    % is harmony?
                                    peaks(peaksIdx,6)=1;
                                    % average amp
                                    %peaks(peaksIdx,8)=average2;
                                    
                                    % peak number in window
                                    peaks(peaksIdx,8)=j;
                                    peaksIdx = peaksIdx +1;
                                    break;
                                end
                            else
                            end
                        end
                        if( octave == 0)
                            'no octave';
                            %start time
                            peaks(peaksIdx,1)=currTime + peaksLoc(j,1)/Fs;
                            %end time
                            peaks(peaksIdx,2)=currTime + peaksLoc(j,2)/Fs;
                            %segment number:
                            peaks(peaksIdx,3)=figureNum;
                            %freq
                            peaks(peaksIdx,4)=newFreq;
                            %ampli
                            peaks(peaksIdx,5)=spectrum(i);
                            % is harmony?
                            peaks(peaksIdx,6)=0;
                            % peak number in window
                            peaks(peaksIdx,8)=j;
                            
                            peaksIdx = peaksIdx +1;
                            if (peaksIdx==7)
                                break
                            end
                        end
                    else continue;
                    end
                end
            end
        else
            break
        end
    end
    
    %   ==========================  Frequency Plot: ========================
    %     figure(10*figureNum+j);
    %     figure;
    if(plotOn == 1)
        
        subplot(1,2,2);
        plot(linspace(0,Fs/2,length(spectrum)),spectrum);
        axis([0 1500 0 1]);
        hold on
        a= size(peaks,1);
        for m = 1 : size(peaks,1)
            if(peaks(m,4)<1500)
                plot(peaks(m,4), peaks(m,5), 'r.');
            end
        end
        hold off
    end
    % look for the desired notes in the notes we found:
    for l=1:length(desiredNotes)
        found = 0;
        if(desiredNotes(l)>0)
            for m=1:size(peaks)
                % if we found the desired note:
                if (desiredNotes(l)+4>peaks(m,4) &&desiredNotes(l)-4<peaks(m,4))
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
            %continuous note?:
            results1(resultIdx,8)=continuousNote;
            % increment index:
            resultIdx=resultIdx+1;
        else
            break
        end
        
    end
    peaks = zeros(6,8);
    peaksIdx = 1;
    prevDesiredNotes = desiredNotes;
    currTime = yEnd/Fs;
    figureNum = figureNum + 1;
    
    ' end';
end
's';
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
'end';

end





