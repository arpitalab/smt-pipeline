function [jumps,vanHove,bins] = calc_vanHove(tracks,lagtime,filter_flag)
% calc_vanHove calculates van Hove correlation at a lag lagtime
% van Hove correlation is defined as G(r,\tau) =
% <delta(r-\vert(r(t+t_0)-r(t)\vert)> where the average is over the initial
% positions of the particle in the trajectory and \delta is the Dirac delta
% function
%
% Input: tracks: cell array with (x,y) positions of tracked particles
%        lagtime: index of desired lag
%        filter_flag : filter immobile_tracks
%        
% Output: 
%        jumps: array of displacements
%        vanHove: vanHove correlation
%        bins: bin centers of the distribution
vH = zeros(1,128);
thresh = 10e-3;
edges=linspace(thresh,(lagtime)*0.4,129); % set up edges of bins to be 1 nm to 416nm*lagtime
jumps=[];
if(nargin<3)
    filter_flag = 0;
end
if(filter_flag)
    immobile_tracks = filter_immobile_tracks(tracks,0.0019,lagtime);
    %sprintf('filtered out %d immobile tracks',length(immobile_tracks))
else
    immobile_tracks = [];
end
jmp_count=0;
for ii=1:length(tracks{1})
    if(length(tracks{1}{ii})>lagtime && isempty(find(immobile_tracks==ii)))
        x = tracks{1}{ii}(:,1); %in um
        y = tracks{1}{ii}(:,2); %in um
        r = sqrt((x(lagtime+1:end)-x(1:end-lagtime)).^2+(y(lagtime+1:end)-y(1:end-lagtime)).^2);
        jumps = [jumps;r];
        if(length(find(r>thresh))>0)
            jmp_count = jmp_count+length(find(r>thresh));
            [vanHove,bins] = histcounts(r(find(r>thresh)),edges); % bin the jumps for each track (discard all jumps that are close to localization precision)
            vH = vH + vanHove/length(find(r>thresh)); % accumulate the histogram after averaging for the number of initial positions
        end
    end
end

bins = 0.5*(bins(2:end)+bins(1:end-1)); % bin centers
vH = vH./(2*pi*bins); % for radial van Hove, need to divide by 2\pi r  ;
vanHove = vH/(trapz(bins,2*pi*bins.*vH)); % convert to PDF by normalizing
%sprintf('%d',jmp_count)

end

