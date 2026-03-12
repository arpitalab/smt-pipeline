function [tracklen,x,f,flo,fup,jump_hist] = find_bound_particles_2(tracks,rmin,rmax)
    unbound={};
    long_jmp={};
    jump_hist=[];
    delta = 7;
    for itrack=1:length(tracks{1})
        xdum=tracks{1}{itrack}(:,1:2);
        if(length(xdum)>delta) % restrict to at least 1.4 sec tracks (e.g. pEM)
            dx=diff(xdum(:,1));
            dy=diff(xdum(:,2));
            dr=sqrt(dx.^2+dy.^2);
            jump_hist = [jump_hist;dr];
            long_jmp{itrack} = find(dr>rmin);
            % need to now keep all track segments that do jump less than rmin
            for iseg=1:length(dx)-delta
                dx=(xdum(iseg+delta,1)-xdum(iseg,1));
                dy=(xdum(iseg+delta,2)-xdum(iseg,2));
                if(sqrt(dx.^2+dy.^2)>rmax)
                    unbound{itrack}(iseg)=1;
                else
                    unbound{itrack}(iseg)=0;
                end
            end
        end
    end
       
    % find all tracks with no large jumps
    good_tracks=cellfun(@isempty,long_jmp);
    % find all tracks with no breaks
    breaks=cellfun(@nnz,unbound);
    if(length(good_tracks) > length(breaks))
        indices=find(good_tracks(1:length(breaks))==1 & breaks==0);
    else
        indices=find(good_tracks==1 & breaks(1:length(good_tracks))==0);
    end
    % get the intersection of these (for now) -- otherwise we have to get the
    % length of the intervening segments as well
    %indices=find(good_tracks==1 & breaks==0);
    bound_tracks={};
    for ii=1:length(indices)
        bound_tracks{ii} = tracks{1}{indices(ii)};
    end
    tracklen=cellfun('size',bound_tracks,1);
    [f,x,flo,fup]=ecdf(tracklen*0.2,'alpha',0.05,'Function','survivor');
    if(x(1)==x(2)) x(1)=x(1)-0.2;
end