% linearSubdivision.m
%
% Copyright (C) 2014 < Olivier DAZEL <olivier.dazel@univ-lemans.fr> >
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


function [newVertices, newFaces] =  linearSubdivision(vertices, faces)
% Linear subdivision for triangle meshes
%
%  Dimensions:
%    vertices: 3xnVertices
%    faces:    3xnFaces
%  
%  Author: Jesus Mena

	global edgeVertex;
    global newIndexOfVertices;
	newFaces = [];
	newVertices = vertices;

	nVertices = size(vertices,2);
	nFaces    = size(faces,2);
	edgeVertex= zeros(nVertices, nVertices);
	newIndexOfVertices = nVertices;

    % ------------------------------------------------------------------------ %
	% create a matrix of edge-vertices and a new triangulation (newFaces).
    % 
    % * edgeVertex(x,y): index of the new vertex between (x,y)
    %
    %  0riginal vertices: va, vb, vc.
    %  New vertices: vp, vq, vr.
    %
    %      vb                   vb             
    %     / \                  /  \ 
    %    /   \                vp--vq
    %   /     \              / \  / \
    % va ----- vc   ->     va-- vr --vc 
	%
    
	for i=1:nFaces
		[vaIndex, vbIndex, vcIndex] = deal(faces(1,i), faces(2,i), faces(3,i));
		
		vpIndex = addEdgeVertex(vaIndex, vbIndex);
		vqIndex = addEdgeVertex(vbIndex, vcIndex);
		vrIndex = addEdgeVertex(vaIndex, vcIndex);
		
		fourFaces = [vaIndex,vpIndex,vrIndex; vpIndex,vbIndex,vqIndex; vrIndex,vqIndex,vcIndex; vrIndex,vpIndex,vqIndex]';
		newFaces  = [newFaces, fourFaces]; 
    end;
    	
    % ------------------------------------------------------------------------ %
	% positions of the new vertices
	for v1=1:nVertices-1
		for v2=v1:nVertices
			vNIndex = edgeVertex(v1,v2);
            if (vNIndex~=0)
 				newVertices(:,vNIndex) = 1/2*(vertices(:,v1)+vertices(:,v2));
            end;
        end;
    end;
 	
end

% ---------------------------------------------------------------------------- %
function vNIndex = addEdgeVertex(v1Index, v2Index)
	global edgeVertex;
	global newIndexOfVertices;

	if (v1Index>v2Index) % setting: v1 <= v2
		vTmp = v1Index;
		v1Index = v2Index;
		v2Index = vTmp;
	end;
	
	if (edgeVertex(v1Index, v2Index)==0)  % new vertex
		newIndexOfVertices = newIndexOfVertices+1;
		edgeVertex(v1Index, v2Index) = newIndexOfVertices;
	end;

	vNIndex = edgeVertex(v1Index, v2Index);

    return;
end