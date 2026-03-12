function [immobile_tracks] = filter_immobile_tracks(tracks,threshold,lagtime)
% immobile_tracks = filter_immobile_tracks(tracks,threshold,lagtime)
% Helper routine identify immobile tracks for H2B before calculating the van Hove
% correlation
% Input: 
%       tracks - cell array of tracks
%       threshold - threshold for msd below which tracks are called
%       immobile (typically 2-4 x localization precision^2
%       lagtime - lag time for calculating immobile tracks
% Output:
%       immobile_tracks - index of all immobile tracks

% now sort tracks based on P(M)
k=1;
immobile_tracks=[];
for ii=1:length(tracks{1})
    if(length(tracks{1}{ii})>7)
        x=tracks{1}{ii}(:,1);
        y=tracks{1}{ii}(:,2);
        rsq = ((x(lagtime+1:end)-x(1:end-lagtime)).^2+(y(lagtime+1:end)-y(1:end-lagtime)).^2);
        mrsq=mean(rsq);
        if(mrsq < threshold)
           immobile_tracks(k) = ii;
           k=k+1;
        end 
    end
end
end

