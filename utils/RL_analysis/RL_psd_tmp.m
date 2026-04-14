function [M,Gs,P1norm,psf] = RL_psd_tmp(bins,vanHove,Niter,sigma_loc)
% Richardson-Lucy algorithm to extract P(D) from G(r,t)
% see Ashwin et al. PNAS
% Given a van Hove correlation function, this routine uses the RL algorithm
% to calculate the underlying distribution of MSDs that gave rise to the vH
% function, with the assumption of Gaussian distributed jumps for a given M
% $G(r,t) = \int dM P(M) A exp(-r^2/M)$
% Input:
%       bins: points at which vanHove correlation is calculated
%       vanHove: vanHove distribution
%       Niter: number of iterations for RL algorithm
%       sigma_loc: (optional) localization precision (µm). When provided,
%                  the PSF kernel is corrected for localization noise:
%                  the observed MSD = true MSD + 4*sigma_loc^2, so the
%                  kernel uses M_eff = M + 4*sigma_loc^2 and the recovered
%                  P(M) represents the TRUE MSD distribution.
%                  Default: 0 (no correction).
% Output:
%       M: values over which empirical MSDs are estimated (TRUE MSDs)
%       Gs: estimated van Hove correlation
%       P1norm: P(M), probability distribution of MSDs

if nargin < 4 || isempty(sigma_loc)
    sigma_loc = 0;
end
noise_var = 4 * sigma_loc^2;   % localization noise contribution to observed MSD

if(size(bins,1)==1)
    x=bins';
else
    x=bins;
end
if(size(vanHove,1)==1)
    vanHove = vanHove';
end
M = logspace(log10(1e-3),log10(1),100); % allocate the grid for MSD calculation (TRUE MSDs)
lM = length(M);
lx = length(x);
psf = zeros(lx,lM);
for ii=1:lM
    M_eff = M(ii) + noise_var;    % observed MSD = true MSD + localization noise
    psf(:,ii) = exp(-(x).^2/M_eff); %Gaussian "PSF" with noise correction
    psf(:,ii) = psf(:,ii)/(pi*M_eff); % normalize PSF
end

P1=exp(-(M)/1e-3); % initial guess. changing the denominator by 2 OoM doesn't make a difference
P1norm=P1/trapz(M,P1); %normalize P(M) so that integral of P1norm is 1

for iterations=1:Niter
    Gs = trapz(M,P1norm(ones(lx,1),:).*psf,2); %RL first step
    Gsest = vanHove./Gs; %ratio of empirical vH to estimated vH
    convest = trapz(x,2*pi*x(:,ones(lM,1)).*Gsest(:,ones(lM,1)).*psf,1); % blur again with PSF
    P1norm = P1norm.*convest; % compute new estimate of P(M)
    P1norm = P1norm/trapz(M,P1norm); %normalize P(M)
    %P1norm = P1norm/sum(P1norm);
    residual=sum((Gs-trapz(M,P1norm.*psf,2)).^2);
    if(residual<1e-12)
        %sprintf('number of iterations = %d',iterations)
        break 
    end
end
Gs = trapz(M,P1norm(ones(lx,1),:).*psf,2); % estimated van Hove correlation
%sprintf('residual %0.5g',residual)
%sprintf('integrated squared error %f',trapz(bins,2*pi*bins'.*((Gs(1:end)-vanHove).^2)))
end
