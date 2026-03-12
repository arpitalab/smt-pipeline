function [x1,msd,msd_err,msd_dist] = msd_calc(classified_tracks,tracks,max_lag)

numgroups = size(classified_tracks,1);
msd = zeros(max_lag,numgroups);
msd_err = zeros(max_lag,numgroups);
min_len=7;
x1=linspace(0.01,0.01*max_lag,max_lag)';
msd_dist = {};
for igroup=1:numgroups
    e_msd = zeros(max_lag,1);
    k=1;
    track_msd={};
    track_std={};
    num_lags={};
    tmp_msd=zeros(max_lag,1);
    tmp_std=zeros(max_lag,1);
    num_lags_tracks=zeros(max_lag,1);
    e_msd=zeros(max_lag,1);
    e_var = zeros(max_lag,1);
    for itrack = 1:length(classified_tracks{igroup})
        L=length(tracks{1}{classified_tracks{igroup}(itrack)});
        if(L>min_len) % only selecting tracks that are longer than 7 frames
            x = tracks{1}{classified_tracks{igroup}(itrack)}(:,1);
            y = tracks{1}{classified_tracks{igroup}(itrack)}(:,2);
            max_lag_track = min(L-1,max_lag);
            for lagtime=1:max_lag_track
                r = ((x(lagtime+1:end)-x(1:end-lagtime)).^2+(y(lagtime+1:end)-y(1:end-lagtime)).^2);
                tmp_msd(lagtime) = sum(r(:)); % compute all squared displacements for the given track
                tmp_std(lagtime) = std(r(:)).^2;
                tmp_std(isnan(tmp_std)) = 0;
                num_lags_tracks(lagtime) = length(r);
            end
            track_msd{k} = tmp_msd(1:max_lag_track);
            track_std{k} = tmp_std(1:max_lag_track);
            num_lags{k} = num_lags_tracks(1:max_lag_track);
            k=k+1; 
        end
    end
    % now compute ensemble mean msd

    dum=zeros(k-1,1);
    for itrack=1:k-1
        dum(itrack) = track_msd{itrack}(4)/num_lags{itrack}(4);
        for lagtime=1:length(track_msd{itrack})
            e_msd(lagtime)=e_msd(lagtime)+track_msd{itrack}(lagtime);
            e_var(lagtime)=e_var(lagtime)+(num_lags{itrack}(lagtime)-1)*track_std{itrack}(lagtime);
            num_lags_tracks(lagtime) = num_lags_tracks(lagtime)+num_lags{itrack}(lagtime);
        end
    end
    msd_dist{igroup} = dum;
    e_msd = e_msd./num_lags_tracks;
    e_var = sqrt(e_var)./(num_lags_tracks-k-1);
    msd(:,igroup) = e_msd;
    msd_err(:,igroup) = 2.576*e_var; %99% confidence interval
end