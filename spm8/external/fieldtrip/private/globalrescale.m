function [H] = globalrescale(f);

% GLOBALRESCALE creates the homogenous spatial transformation matrix
% for a 7 parameter rigid-body transformation with global rescaling
%
% Use as
%   [H] = globalrescale(f)
%
% The transformation vector f should contain the 
%   x-shift
%   y-shift
%   z-shift
% followed by the
%   pitch (rotation around x-axis)
%   roll  (rotation around y-axis)
%   yaw   (rotation around z-axis)
% followed by the
%   global rescaling factor

% Copyright (C) 2000-2005, Robert Oostenveld
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

% This file is part of FieldTrip, see http://www.ru.nl/neuroimaging/fieldtrip
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id: globalrescale.m 7123 2012-12-06 21:21:38Z roboos $

% compute the homogenous transformation matrix for the translation
T = translate(f([1 2 3]));

% compute the homogenous transformation matrix for the rotation
R = rotate(f([4 5 6]));

% compute the homogenous transformation matrix for the global scaling
S = scale(f([7 7 7]));

% compute the homogenous transformation matrix for the combination
H = T*R*S;

