function f = getLSField2(f, t, lsType, resPhi, normFlag, varargin)
% ==============================================================================
% optional: tile - integer, unit cell periodicity
%           L    - length of unit cell sides (currently cubic only)
% 
% Chan, Y.-C., Ahmed, F., Wang, L., and Chen, W., 2020.
% "METASET: Exploring Shape and Property Spaces for Data-Driven Metamaterials Design.“
% Journal of Mechanical Design. March 2021; 143(3): 031707.
% 
% Author: Yu-Chin Chan (ychan@u.northwestern.edu)
% Last updated: 4/28/2020
%               9/10/2020 - added normalization option
% ==============================================================================
if nargin == 6
    tile = varargin{1}; L = 1; %L = [1,1,1];
elseif nargin == 7
    [tile, L] = varargin{:}; 
else
    tile = [1 1 1]; L = 1; %L = [1,1,1];
end
X = linspace(0,2*pi*tile(1)/L,tile(1)*resPhi); 
Y = linspace(0,2*pi*tile(2)/L,tile(2)*resPhi);
Z = linspace(0,2*pi*tile(3)/L,tile(3)*resPhi);
[X,Y,Z] = meshgrid(X,Y,Z);
f = f(X,Y,Z);
if normFlag, f = f./max(abs(f(:))); end
switch lsType
    case 'tm', f = f - t;
    case 'tp', f = t - f;
    case 'tt', f = f.^2 - t.^2;
end
end