function [M,Gs,P1norm,psf] = RL_psd_tmp(bins,vanHove,Niter)
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
% Output:
%       M: values over which empirical MSDs are estimated
%       Gs: estimated van Hove correlation
%       P1norm: P(M), probability distribution of MSDs
% A test case is the empirical distribution of 2 Gaussians 
% Memp=logspace(log10(1e-3),log10(2e-1),200);
% PM=30*exp(-(M-5e-3).^2/(2*1e-6))+10*exp(-(M-5e-2).^2/(2*1e-4));
% which can be used to forward calculate G(r,t) using the point-spread
% function below. Then RL can be applied to recover this "unknown" G(r,t) as
% below.

if(size(bins,1)==1)
    x=bins';
else
    x=bins;
end
if(size(vanHove,1)==1)
    vanHove = vanHove';
end
M = logspace(log10(1e-3),log10(1),100); % allocate the grid for MSD calculation
lM = length(M);
lx = length(x);
psf = zeros(lx,lM);
for ii=1:lM
    psf(:,ii) = exp(-(x).^2/(M(ii))); %Gaussian "PSF"
    psf(:,ii) = psf(:,ii)/(pi*M(ii)); % normalize PSF
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
