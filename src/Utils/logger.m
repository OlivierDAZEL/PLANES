% logger.m
%
% Copyright (C) 2016 Mathieu GABORIT <gaborit@kth.se>
%
% This file is part of PLANES.
%
% PLANES (Porous LAum NumErical Simulator) is a software to compute the
% vibroacoustic response of sound packages containing coupled
% acoustic/elastic/porous substructures. It is mainly based on the
% Finite-Element Method and some numerical methods developped at
% LAUM (http://laum.univ-lemans.fr).
%
% You can download the latest version of PLANES at
% https://github.com/OlivierDAZEL/PLANES
% or find more details on Olivier's webpage
% http://perso.univ-lemans.fr/~odazel/
%
% For any questions or if you want to
% contribute to this project, contact Olivier.
%
% PLANES is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program. If not, see <http://www.gnu.org/licenses/>.
%%

function logger(verbosity, level, section, msg, wait)
	% logging function for verbous runs
	%
	% verbosity -- current verbosity level
	% level -- verbosity level required to show message
	% section -- section of the programm trigging the logger
	% msg -- message

	if level==0
		prefix = '=>';
		out = stdout;
	else
		prefix = '[log]';
		out = stderr;
	end
	if verbosity>=level
		fprintf(out, '%s %s: %s\n', prefix, section, msg);
	end
end
