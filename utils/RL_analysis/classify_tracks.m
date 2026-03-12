function [classified_tracks] = classify_tracks(tracks,P1norm,M,lagtime,iplot,title_str)
%UNTITLED5 Summary of this function goes here
%   Detailed explanation goes here

% now sort tracks based on P(M)
jdx=islocalmin(P1norm); % find all local minima in the P(M) function. these will serve as classification boundaries
idum = find(jdx);
nminima = length(idum);
classified_tracks=cell(nminima+1,1);
if(nminima>0) 
    for ii=1:length(tracks{1})
        if(length(tracks{1}{ii})>5)
            x=tracks{1}{ii}(:,1);
            y=tracks{1}{ii}(:,2);
            rsq = ((x(lagtime+1:end)-x(1:end-lagtime)).^2+(y(lagtime+1:end)-y(1:end-lagtime)).^2);
            mrsq=mean(rsq);
            m=1;
            if(mrsq < M(idum(1)))
                classified_tracks{m} = [classified_tracks{m};ii]; % accumulate all tracks in this range
            end
            if(nminima>1)
                for im=1:nminima-1
                    if(mrsq > M(idum(im)) && mrsq < M(idum(im+1)))
                        classified_tracks{m+im} = [classified_tracks{m+im};ii]; % accumulate all tracks in this range  
                    end
                end
            end
            m = nminima+1;
            if(mrsq > M(idum(nminima)))
                classified_tracks{m} = [classified_tracks{m};ii]; % accumulate all tracks in this range
            end
        end  
    end
end
if(iplot)
    figure;
    sgtitle(title_str);
    for ii=1:nminima+1
        subplot(1,nminima+1,ii);
        for jj=1:min(500,length(classified_tracks{ii}))
            x=tracks{1}{classified_tracks{ii}(jj)}(:,1);
            y=tracks{1}{classified_tracks{ii}(jj)}(:,2);
            plot(x,y);
            hold on;
        end
        axis equal;axis([5 10 5 10]);
        title(num2str(length(classified_tracks{ii})/length(tracks{1})));
    end
end
end

