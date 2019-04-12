function [psi,X] = cmorwavf(LB,UB,N,cmrl)
%CMORWAVF Complex Morlet wavelet.
%   [PSI,X] = CMORWAVF(LB,UB,N,FB,FC) returns values of
%   the complex Morlet wavelet defined by a positive bandwidth
%   parameter FB, a wavelet center frequency FC, and the expression
%   PSI(X) = ((pi*FB)^(-0.5))*exp(2*i*pi*FC*X)*exp(-(X^2)/FB)
%   on an N point regular grid in the interval [LB,UB].
%
%   Output arguments are the wavelet function PSI
%   computed on the grid X.
%
%   See also WAVEINFO.

%   M. Misiti, Y. Misiti, G. Oppenheim, J.M. Poggi 09-Jun-99.
%   Last Revision: 14-May-2003.
%   Copyright 1995-2004 The MathWorks, Inc.
%   $Revision: 1.8.4.2 $  $Date: 2004/03/15 22:39:54 $


Fc = 1; Fb = 1;

% Compute values of the Complex Morlet wavelet.
X = linspace(LB,UB,N);  % wavelet support.
psi = ((pi*Fb)^(-0.5))*exp(2*i*pi*Fc*X).*exp(-(X.*X)/Fb);

%VOOR GEBRUIK:
%wavemngr('del', 'cmrl')
%wavemngr('add', 'CMorlet', 'cmrl', 5, '', 'cmorwavf', [-3, 3])