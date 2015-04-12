% FINDFREQUENCIES Takes a recording of a piano and returns the frequencies
% of the notes that the piano played.
%
% Inputs:
%        song = Matlab vector of the wav file recording
%          fs = Sampling frequency of the wav file
% notechanges = Vector of song indices where the notes being played change
%               (from timeout.m)
%
% Output:
% frequencies = Vector of the frequencies being played (length is one
%               shorter than notechanges)

% Scott Steger
% ELEC 301
% 11 December 2006

function frequencies = findfrequencies(song,fs,notechanges)

    % Normalize song
    song = song/max(abs(song));
    
    % Initialize output in case no notes are played at all
    frequencies = zeros(length(notechanges)-1,1);

    % Repeat for each different segment of notes as identified by timeout.m
    for(numSegment = 1:length(notechanges)-1)
        
        % Extract only the segments of the song when the notes don't change
        songSegment = song(notechanges(numSegment):notechanges(numSegment+1));
        
        % Take the Fourier transform of the segment and limit it to
        % non-redundant frequencies
        spectrum = fft(songSegment,8*length(songSegment));
        spectrum = abs(spectrum(1:ceil(length(spectrum)/2)));

        numNote = 1;
        
        % Assume a note is being played if any spectrum peak is above 500
        if(max(spectrum) > 200)
            
            %filter = gaussfilt(128);
            spectrum = spectrum .^ 2;
            %spectrum = conv(spectrum,filter);
            spectrum = spectrum/max(spectrum);
            
            % Uncomment to plot the FFT of each segment.
            figure(numSegment)  
            plot(linspace(0,fs/2,length(spectrum)),spectrum)
            axis([0 1500 0 1])
            
            % Find the peak of each harmonic that is above a certain threshold
            for(i = 2:length(spectrum)-1)
                
                if(spectrum(i) > .1 & spectrum(i) > spectrum(i-1) & spectrum(i) > spectrum(i+1))
                    
                    % Do not record dual-peaks by only incrementing the counter if the newest detected peak is significantly different than the previous
                    if(abs(i*fs/(length(spectrum)*2) - frequencies(numSegment,numNote)) > 10)
                        
                        newFrequency = i*fs/(length(spectrum)*2);
                        spectrum(i);
                        
                        % Check if the peak is a possible strong harmonic and remove it if it is
                        if(numNote >= 2)
                            possibleHarmonics = 0;
                            
                            % Iterate over each previously determined frequency
                            for(j = 2:numNote)
                                
                                % The note is a possible harmonic if it's frequency is within 3% of the perfect harmonic
                                if(mod(newFrequency,frequencies(numSegment,j)) < .03*newFrequency | ...
                                    frequencies(numSegment,j) - mod(newFrequency,frequencies(numSegment,j))  < .03*newFrequency)
                                    
                                    % If the magnitude of the peak is strong enough, declare that it is a different note instead of a loud harmonic
                                    if(spectrum(i) < .3*frequencies(numSegment,j)*2*length(spectrum)/fs)
                                        possibleHarmonics = 1;
                                    end
                                end
                            end
                            if(~possibleHarmonics)
                                numNote = numNote + 1;
                                frequencies(numSegment,numNote) = newFrequency;
                            end
                        else
                            numNote = numNote + 1;
                            frequencies(numSegment,numNote) = newFrequency;
                        end
                    end
                end
            end
        end
    end
    
    % If any notes have been played, remove the default output
    [rows,columns] = size(frequencies);
    if(columns >= 2)
        frequencies = frequencies(:,2:columns);
    end
end

% ------------- % This is code to make the edge detecting filter % ----%
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