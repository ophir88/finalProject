% CONVERTSONG is a wrapper for all of the functions needed to convert a wav
% file to the notes that are being played
%
% Inputs:
%  filename = String containing the location of the wav file to be
%             converted.
%
% Outputs:
%     music = Letter names of the notes that were played during each
%             interval, and the length of each interval.
%     tempo = The rate at which the song was played.

% Scott Steger
% ELEC 301
% 11 December 2006

function [music,tempo] = convertsong(filename)

    [song,fs] = audioread(filename);
    points = timeout(song,fs);
    frequencies = findfrequencies(song,fs,points);
    notes = frequencyToNote(frequencies);
    [music,tempo] = noteDuration(points,notes,fs);

end