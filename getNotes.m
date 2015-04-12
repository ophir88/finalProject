function[]  = getNotes(Notes,noteIdx,NoteCount, currTime )


    %     currTime*Fs
    %     L
    
    % -------------------------------------------- new loop:
    
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
    
    
    
    if((currTime+windowSize/Fs) > Notes(noteIdx,5)+0.001)
        yEnd = Notes(noteIdx,5)*Fs;
    else
        yEnd = currTime*Fs+windowSize;
    end
    if (yEnd > L)
        yEnd = L;
    end
    desiredNotes = tempDesiredNotes;
    start = currTime;
    dend = yEnd/Fs;
    'a';
    % -------------------------------------------- end new loop
    
    
    
    
    
    %     desiredNotes = zeros(6,1);
    
    %     theEnd=0;
    ContinousNote =0;
    %     if (k+windowSize>L)
    %         partSignal = mySignal(k:L)
    %         theEnd = 1;
    %     else
    %         partSignal = mySignal(k:k+windowSize);
    %     end
    
    
    
    % ****************** get desired notes: ******************
    %set note we look for:
    %     j=1;
    %     for i=notesIdx:length(Notes)
    %
    %         if j<7
    %             if (round(Notes(notesIdx,5)/windowSizeTime)==notesSegIdx)
    %                 desiredNotes(j)=Notes(notesIdx,3);
    %                 notesIdx=notesIdx+1;
    %                 j=j+1;
    %
    %             else
    %                 j=7;
    %                 break;
    %             end
    %         else
    %             break
    %         end
    %
    %
    %     end
    %     notesSegIdx=notesSegIdx+1;
    
    
    
    % ****************** find envelope ******************
    
    % Normalize the sound sample:
    's';
    if(currTime==0)
        currTime = 0.0001;
    end
    yStart = round(currTime*Fs);
    sound(mySignal(yStart:yEnd),Fs)
    % before absolute value, keep original for fourier transform:
    yFourier = mySignal(yStart:yEnd);
    y = abs(yFourier);
    ylength = length(y)/Fs;
    maxenv = y;
    peakwindowsize = 512;
    % get envelop vector:
    average = 0;
    's';
    %     maxenvaa=length(maxenv)
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
    
    average = average/length(y);
    % find amp peaks:
    peaksLoc = zeros(4,2);
    inPeak = 0;
    peakStartLoc = 0;
    %     peakEnd = 0;
    peakEndLoc= 0;
    thresh = average;
    peakInd =1;
    continuousNote = 0;
    for i=1:length(maxenv)
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
                if(peakEndLoc-peakStartLoc>windowSize/8)
                    
                    peaksLoc(peakInd,1)=peakStartLoc;
                    peaksLoc(peakInd,2)=peakEndLoc;
                    peakInd = peakInd+1;
                end
            end
        end
    end
    if(inPeak == 1 && peakStartLoc<length(y)*7/8)
        peaksLoc(peakInd,1)=peakStartLoc;
        peaksLoc(peakInd,2)=length(y);
        continuousNote = 1;
    end
    
    
    % %
    figure(figureNum+101)
    plot (maxenv, 'm')
    %
    for i=1:peakInd
        if (peaksLoc(i,2)-peaksLoc(i,1)>0)
            hold on
            plot([peaksLoc(i,1) peaksLoc(i,2)],[thresh thresh],'--or')
            
        else
            break
        end
    end
    
    %     check envelope of half note + two quarter notes at the same time.
    
    continuousNote;
    's';
    
    
    
    % loop over peak locations, looking for frequencies:
    
    for j=1:peakInd
        
        if (peaksLoc(j,2)-peaksLoc(j,1)>0)
            %note length in quarters:
            %             noteLength = round((peaksLoc(j,2)-peaksLoc(j,1))/(windowSize/2));
            's';
            %-*-*-*----------------**-*-*-*-* fourier etc.*-**-*--**--**--*-*+9++**+**++*+*
            
            NFFT = 2^nextpow2(peaksLoc(j,2)-peaksLoc(j,1)); % Next power of 2 from length of y
            
            spectrum = fft(yFourier(peaksLoc(j,1):peaksLoc(j,2)),NFFT);
            spectrum = abs(spectrum(1:ceil(length(spectrum)/2)));
            
            % max(spectrum);
            filter = gaussfilt(128);
            spectrum = spectrum .^ 2;
            % regular fft:
            spectrum = spectrum/max(spectrum);
            tranCoef = Fs/(length(spectrum)*2);
            
            %             plot signals:
            figure(100*figureNum+j)
            %             %                 figure(1)
            %             %
            plot(linspace(0,Fs/2,length(spectrum)),spectrum)
            axis([0 1500 0 1])
            %     figure(round(k/windowSize)+11)
            %     plot(linspace(0,Fs/2,length(spectrum2)),spectrum2)
            %     axis([0 1500 0 1])
            
            % ------------------- find peaks: ----------------
            
            % look only at relevant frequencies (note frequency spectrum)
            'loop';
            for i=round(27*length(spectrum)/Fs):round(4180*length(spectrum)/Fs)
                
                %variable for close peaks:
                minorPeak = 0;
                % max point:
                if(spectrum(i)>threshold & spectrum(i)>spectrum(i-1)& spectrum(i)>spectrum(i+1))
                    
                    if (peaksIdx>6)
                        break;
                    end
                    % make sure peaks are not too close
                    
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
                                    
                                    %                             'added'
                                    break;
                                end
                                
                            else
                            end
                            
                            
                        end
                        if( octave == 0)
                            
                            'no octave';
                            
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
                        %                 end
                    end
                    
                end
            end
            % ------------------- end find peaks -------------
            
            %-*-*-*----------------**-*-*-*-* fourier etc.*-**-*--**--**--*-*+9++**+**++*+*
            
            
            
        else
            break
        end
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
    pastTime = currTime;
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




