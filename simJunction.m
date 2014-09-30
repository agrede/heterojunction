function [eta,jn] = simJunction(Stack,T,Param,PC,psisRng,neta,approxC,approxV)
% SIMJUNCTION simulates electrostatics at equalibrium by solving the carrier
%   concentrations and poisson equation
%
%   [ETA, JN] = SIMJUNCTION(STACK,T,PROP,PC)
%       STACK   = gate stack struct (must use SI units)
%       T       = temperature in K
%       PROP    = material property parameters (band params, dielectrics, etc.)
%       PC      = physical constants (epsilon0, boltzmann constant, etc.)
%   [ETA, JN] = SIMJUNCTION(STACK,T,PROP,PC,PSISRNG)
%       PSISRNG = range of simulaiton in eV [|below val band|,|above cond band|]
%                       (default [0.5 0.5])
%   [ETA, JN] = SIMJUNCTION(STACK,T,PROP,PC,PSISRNG,NETA)
%       NETA    = number of points across range to calculate (default 5001)
%   [ETA, JN] = SIMJUNCTION(STACK,T,PROP,PC,PSISRNG,NETA,APPROXC,APPROXV)
%       APPROXC = approximations to use for [Gamma;X;L] valley in electron conc
%                       (default [3;2;2])
%       APPROXV = approximations to use for [hh;lh;so] in hole conc
%                       (default [2;2;2])
%       Approximations: (Note: some settings could give wildly inaccurate values)
%               0: Ignore this component
%               1: Maxwell-Boltzmann with parabolic band
%               2: Fermi-Dirac with parabolic band
%               3: Fermi-Dirac with non-parabolic band
%
% Copyright (C) 2014 Alex J. Grede
% GPL v3, See LICENSE.txt for details
% This function is part of heterojunction (https://github.com/agrede/heterojunction)

  if (nargin < 6)
    psisRng = [0.5 0.5];
  endif
  if (nargin < 7)
    neta = 5001;
  endif
  if (nargin < 8)
    approxC = [3;2;2];
    approxV = [2;2;2];
  endif

  kT = PC.kB.*T;

  jn = cell;

  Eg = zeros(length(Stack),1);
  VBO = zeros(length(Stack),1);
  for k = 1:length(Stack)
    jn{k} = semiProps(Stack{k},T,Param,PC);
    Eg(k) = min(jn{k}.Eg);
    jn{k}.Egm = Eg(k);
    jn{k}.etaV = -Eg(k)./kT;
    VBO(k) = jn{k}.VBO;
  endfor

  % reference junction
  [Egr, kr] = max(Eg);

  eta = linspace(-(psisRng(1).*PC.e+Egr),psisRng(2).*PC.e,neta)'./kT;
  detac = (Eg-Egr+VBO-VBO(kr))./kT;

  for k = 1:length(Stack)
    jn{k}.detac = detac(k);
    jn{k}.NI = impurities(eta-detac(k), jn{k}.etaV, Stack{k}.impurities, T,...
                          jn{k}.Impurities, PC);
    jn{k}.n = carConc(eta-detac(k), jn{k}.me, (Eg(k)-jn{k}.Eg)./kT, jn{k}.Eg,...
                      approxC, T, PC);
    jn{k}.p = carConc(-(eta-detac(k)), jn{k}.mh, ...
                      ([0;0;-jn{k}.delta_so]-Eg(k))./kT, ...
                      ([0;0;jn{k}.delta_so])+Eg(k), approxV, T, PC);
    jn{k}.rho = PC.e.*sum([jn{k}.p -jn{k}.n jn{k}.NI], 2);
    [jn{k}.D, jn{k}.C, jn{k}.psi] = poissonSolution((eta-detac(k)),jn{k}.rho,...
                                                    T, jn{k}.kappas, PC);
    jn{k}.phin = -interp1(jn{k}.psi, (eta-detac(k)), 0, 'spline').*kT./PC.e;
    jn{k}.phip = Eg(k)./PC.e-jn{k}.phin;
    [jn{k}.x, jn{k}.ind] = yByEta(jn{k}.psi, jn{k}.D, jn{k}.kappas, PC);

    if k>1
      jn{k-1}.eta0 = interp1(sum([jn{k-1}.D jn{k}.D],2),eta,0,'spline');
      [tmp, jn{k-1}.ind0] = min(abs(eta-jn{k-1}.eta0));
    endif
  endfor
endfunction
