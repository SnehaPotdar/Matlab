function w = decompressShape(shape)
%decompressShape Decompress a gradient or pulse shape.
%   w=decompressShape(shape) Decompress the shape compressed with a run-length
%   compression scheme on the derivative. The given shape is structure with
%   the following fields:
%     num_samples - the number of samples in the uncompressed waveform
%     data - containing the compressed waveform
%
%   See also compressShape
        
T = 1e-8;
dataPack = shape.data;
dataPack = dataPack - mod(dataPack,T);       %To check for meaningful differences between amplitudes
numSamples=shape.num_samples;

w= zeros(1, numSamples) ;                 % eine Matrix wird erstellt, in die später der unkonvertierte Datensatz eingetragen werden kann
                                                % Dimension der Matrix: (1,Länge des gesamten Datensatzes)
                                       
%%%% Dekomprimierung

countPack= 1;                                               % Zähler 1: inkrementiert komprimierten Datensatz
countUnpack= 1;                                             % Zähler 2: inkrementiert unkomprimierten Datensatz      


while countPack < length(dataPack)                          % while-Schleife wird solange abgearbeitet, bis Zähler 1  
%          countPack= round(countPack);                                                  % nachfolgendem ~ verglichen, wenn ungleich, dann ...
%          countUnpack = round(countUnpack);                                             % am Ende des dataPack angekommen ist
    
    if dataPack(countPack) ~= dataPack(countPack + 1)       % aktueller Zahlenwert im komprimierten Datensatz wird mit 
%          countPack= round(countPack);                                                  % nachfolgendem ~ verglichen, wenn ungleich, dann ...
%          countUnpack = round(countUnpack);
        

%     if (abs(dataPack(countPack) - dataPack(countPack + 1)) > T)
        w(countUnpack)= dataPack((countPack));                % ...wird die aktueller Zahl an die entsprechende Stelle im
                                                            % unkomprimierten Datensatz geschrieben
        countUnpack= countUnpack + 1;                       % abschließend Erhöhung im komprimierten und
        countPack= countPack + 1;                           % im unkomprimierten Datensatz um 1
    
    else
                                                            % wenn die beiden Werte hingegen gleich sind, dann ...
        rep= (dataPack(countPack + 2) + 2);                   % rep = Variable, die angibt, wie oft der jeweilige Wert 
        if(mod(rep,1)~=0)
             w(countUnpack)= dataPack((countPack));                % ...wird die aktueller Zahl an die entsprechende Stelle im
            countUnpack= countUnpack + 1;                       % abschließend Erhöhung im komprimierten und
            countPack= countPack + 1;  
        else
                                                                % geschrieben werden muss
        w(countUnpack:(countUnpack + rep - 1))= dataPack(countPack);
                                                            % 
        countPack= countPack + 3;
        countUnpack= countUnpack + rep;
        end
      
    end
end
if countPack==length(dataPack)                              % Last sample is not a repeat, so include it directly
    w(countUnpack)=dataPack(countPack);
end

w = cumsum(w);
w=w(:);