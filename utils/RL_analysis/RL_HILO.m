function [M,bins,computed_quantities] = RL_HILO(tracks,title_str,lagtime)
% Driver function for the Richardson Lucy algorithm
% Inputs:
%       tracks: cell array of tracks (in Kaustubh's format)
%       title_str: string for labeling the plots (e.g. H2B)
%       lagtime: lag to calculate RL analysis 
%
% Outputs:
%       M: vector of MSD values 
%       bins: vector of displacement values for van Hove correlation
%       computed_quantities: struct containing P(MSD), vHC, vector of state
%       classification, msds of each track in each state, fits of msd/lag
%       time to power law
%
% Example run
% tracks = summary_table.X(1); % replace 1 with whatever index you wish to
% analyse
% title_str = extractAfter(summary_table.cell_protein{1},['D4','_']); % you
% can change this to any substring etc. This is just used for labeling the
% figures and nothing else.

%title_str = extractAfter(T.cell_protein{index},[cell_string,'_']);
%if(contains(title_str,"_12ms"))
%    title_str = extractBefore(title_str,"_12ms");
%end
sprintf('computing for %s',title_str)
%tracks=T.X(index);
[~,vanHove,bins]=calc_vanHove(tracks,lagtime,1);
[M,Gs,P1norm,~] = RL_psd_tmp(bins,vanHove,50000);
figure;
plot(bins,vanHove,'ko','markersize',8);
hold on;
set(gca,'yscale','log');
plot(bins,Gs,'r','linewidth',2);
set(gca,'linewidth',1,'fontsize',24,'fontweight','bold');
xlabel('r [\mum]');
ylabel('G(r,\tau)');
axis([0 1 1e-4 100]);
axis square;
box off;
title(title_str);
figure;
semilogx(M,P1norm,'k','linewidth',2);
title(title_str);
axis([1e-3 2e-1 0 200]);
xlabel('MSD [\mum^2]');
ylabel('P(MSD)');
set(gca,'linewidth',1,'fontsize',24,'fontweight','bold');
box off;
classified_tracks=classify_tracks(tracks,P1norm,M,lagtime,1,title_str); %classify tracks based on the RL P(MSD) distribution
[~,msd,msd_err] = msd_calc(classified_tracks,tracks,25); % calculate ensemble msd for the different groups
%
% Fit the ensemble msd data for each group (provided there are
% more than some number of tracks) to a power law model and return that
lsqOpts = optimoptions('lsqnonlin', 'MaxFunctionEvaluations', 1000, 'Display', 'none');
x1 = linspace(0.2,25*0.2,25)'; % this needs to be modified to accept the sampling time and maximum lag as input
fits = {};
fitci = {};
for ii=1:size(msd,2)
   if(length(classified_tracks{ii})>200) % if there are sufficient number of tracks
       dat1 = msd(:,ii);
       Fsumsquares = @(x)(log(fun_v(x)) - log(dat1));
       [fitpars, ResNorm, Residual, Flag, Output, Lambda, Jacobian] = lsqnonlin(@(x) my_fun(x,(1:25)',dat1(:)),[1e-4,1e-4,0.5],[1e-4,1e-6,0.1],[5e-1,10e-3,1],lsqOpts);
       %[fitpars, ResNorm, Residual, Flag, Output, Lambda, Jacobian] = lsqnonlin(Fsumsquares,[0.007,0.4,0.001],[1e-4,0.1,1e-3],[2e-1,0.7,10e-3],lsqOpts);
       CI = nlparci(fitpars,Residual,'jacobian',Jacobian);
       fits{ii}=fitpars;
       fitci{ii}=CI;
   end
end

% Returning a structure with different fields so
% one does not need to compute everything all over again.
computed_quantities.ident = title_str; %identifier
computed_quantities.lagtime = lagtime; %lagtime for computation
computed_quantities.x = bins; %vector of displacements for vHC
computed_quantities.vanHove = vanHove; % van Hove correlation
computed_quantities.P1norm = P1norm; % P(M)
computed_quantities.Gs = Gs; % estimated vHC
computed_quantities.classified_tracks = classified_tracks; %cell array of classified tracks
computed_quantities.msd = msd; % mean squared displacement of all tracks for all groups
computed_quantities.msderr = msd_err; % error in msd
computed_quantities.fits = fits; %msd fits
computed_quantities.fitci = fitci; % confidence interval of fits

end

function F = my_fun(x,lags,e_msd)
  
b = (abs(1+0.05./lags).^(2+x(3))+abs(1-0.05./lags).^(2+x(3))-2)./((0.05./lags).^2);

  mlfit = 2*x(1)/((1+x(3))*(2+x(3)))*(0.2.^x(3)*lags.^x(3).*b-2*0.01.^x(3))+2*x(2);
  
  F = (log(mlfit)-log(e_msd)).^2;
end

