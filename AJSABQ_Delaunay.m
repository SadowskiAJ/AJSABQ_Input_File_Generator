function [TRI,VER,TRI_boundary,VER_boundary] = AJSABQ_Delaunay(BX,BY,BID,IX,IY,IID,tolRND,tolNN,tolBG,order,offsetN,offsetE,scale)
% Constrained Delaunay triangulation algorithm for the 
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 12:47 (previously 01/07/12 - 20:40)

% BSD 3-Clause License
% 
% Copyright (c) 2022, Dr Adam Jan Sadowski of Imperial College London
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions are met:
% 
% 1. Redistributions of source code must retain the above copyright notice, this
%    list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright notice,
%    this list of conditions and the following disclaimer in the documentation
%    and/or other materials provided with the distribution.
% 
% 3. Neither the name of the copyright holder, nor of Imperial College, nor the names of its
%    contributors may be used to endorse or promote products derived from
%    this software without specific prior written permission.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
% AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
% IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
% DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
% FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
% DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
% SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
% CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
% OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
% OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
% 
% THE USE OF THE MATLAB LANGUAGE DOES NOT IMPLY ENDORSEMENT BY MATHWORKS.
% THE USE OF THE MATLAB_MPI TOOLBOX DOES NOT IMPLY ENDORSEMENT BY MIT.

%%%%%%%%%%%%%%%%% Preliminaries and definitions %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     
BOUNDARY = struct('x',[],'y',[],'ID',[]);
% .x - x-coordinate of boundary vertex
% .y - y-coordinate of boundary vertex
% .ID - label of the node corresponding to the vertex for ABAQUS

INTERIOR = struct('x',[],'y',[],'ID',[]);
% .x - x-coordinate of interior vertex
% .y - y-coordinate of interior vertex
% .ID - label of the node corresponding to the vertex for ABAQUS

VER = struct('x',[],'y',[],'boundary',[],'triangle',[],'pos',[],'ID',[]); % Vertex struct
% .x - x-coordinate of vertex
% .y - y-coordinate of vertex
% .boundary - is this a boundary vertex (1 or 0)?
% .triangle - which triangles does the vertex belong to (TRI.pos)?
% .pos - index of vertex within the STRUCT
% .ID - label of the node corresponding to the vertex for ABAQUS

TRI = struct('vertices',[],'edges',[],'pos',[],'ID',[]); % Triangle struct
% .vertices - 3 vertices of the triangle (VER.pos)
% .edges - 3 edges of the triangle (EDGE.pos)
% .pos - index of the triangle within the STRUCT
% .ID - label of the element corresponding to the triangle for ABAQUS

EDGE = struct('vertices',[],'triangle',[],'boundary',[],'pos',[]); % Edge struct
% .vertices - 2 vertices of the edge (VER.pos)
% .triangle - which triangle does the edge belong to (TRI.pos)?
% .boundary - is this a boundary edge (1 or 0)?
% .pos - index of the edge within the struct

% Scaling of input data for improved numerical stability
if scale == 1
    mxX = max(abs(BX)); mxY = max(abs(BY));
    factX = 1.0/mxX; factY = 1.0/mxY;
    BX = BX*factX; IX = IX*factX;
    BY = BY*factY; IY = IY*factY;
    BX = round(BX*10^tolRND)/(10^tolRND); BY = round(BY*10^tolRND)/(10^tolRND);
    IX = round(IX*10^tolRND)/(10^tolRND); IY = round(IY*10^tolRND)/(10^tolRND);
else
    factX = 1.0; factY = 1.0;
end

% Establishing the correct arrays
BOUNDARY.x = BX; BOUNDARY.y = BY; BOUNDARY.ID = BID;
INTERIOR.x = IX; INTERIOR.y = IY; INTERIOR.ID = IID;
VER.x = BOUNDARY.x; 
VER.y = BOUNDARY.y;
VER.boundary = ones(1,length(VER.x));
VER.pos = [1:length(VER.x)]; 
VER.ID = BOUNDARY.ID;
warning off MATLAB:divideByZero

%%%%%%%%%%%%%%%%% Algorithm 1.1 - Protruding point removal - page 11 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if order == 1; disp('Algorithm 1/4 - Boundary triangulation'); end
if order == 2; disp('Algorithm 1/5 - Boundary triangulation'); end

% create an initial triangulation of the boundary
POLYGON.x = VER.x; POLYGON.y = VER.y; POLYGON.boundary = VER.boundary; POLYGON.pos = VER.pos; POLYGON.ID = VER.ID;
EDGE.triangle = zeros(2*length(VER.ID) - 3,2); % Using Eq. 1.3 p. 8
VER.triangle = zeros(length(VER.ID),2*length(INTERIOR.ID) + length(VER.ID) - 2); % Using Eq. 1.2 p. 8
TRI.edges = zeros(length(VER.ID)-2,3); % Using Eq. 1.2 p. 8
% See Oyvind Helle & Morten Daehlen (2006) - Triangulations and Applications - Springer

while length(POLYGON.ID) > 3; choices = []; BIG_PID = []; BIG_Px = []; BIG_Py = []; h = []; A = []; Qk = []; in = []; 
    BIG_PID = [ [length(POLYGON.ID) 1:(length(POLYGON.ID)-1)]' [1:length(POLYGON.ID)]' [2:length(POLYGON.ID) 1]'];
    BIG_Px  = POLYGON.x(BIG_PID); BIG_Py = POLYGON.y(BIG_PID);

    % test 1 - check internal angle (must be < 180 degrees)
    a = [BIG_Px(:,2) - BIG_Px(:,3), BIG_Py(:,2) - BIG_Py(:,3)];
    b = [BIG_Px(:,2) - BIG_Px(:,1), BIG_Py(:,2) - BIG_Py(:,1)];
    ang = real(acos(dot(a,b,2)./([sqrt(a(:,1).^2 + a(:,2).^2)].*[sqrt(b(:,1).^2 + b(:,2).^2)])));

    % test 2 - check if element centroid is inside the polygon
    xC = mean(BIG_Px,2); yC = mean(BIG_Py,2);
    in = EXTERNAL_inpoly([xC yC],[POLYGON.x' POLYGON.y']);

    % test 3 - check that the three points are not collinear
    Col = abs(sum(cross([BIG_Px(:,1) BIG_Px(:,2) BIG_Px(:,3)],[BIG_Py(:,1) BIG_Py(:,2) BIG_Py(:,3)]),2));
    
    % test 4 - check that no more than 3 points of the boundary are inside the triangle
    count = [];
    for S = 1:length(BIG_PID)
        count(end+1) = sum(EXTERNAL_inpoly([BOUNDARY.x' BOUNDARY.y'],[BIG_Px(S,:)' BIG_Py(S,:)']));
    end

    % calculate element quality - Eq. 1.8 - see George & Borouchaki (1998) - Delaunay Triangulation & Meshing - Editions Hermes
    A = 0.5*abs( (BIG_Px(:,2).*BIG_Py(:,1) - BIG_Px(:,1).*BIG_Py(:,2)) + (BIG_Px(:,3).*BIG_Py(:,2) - BIG_Px(:,2).*BIG_Py(:,3)) + (BIG_Px(:,1).*BIG_Py(:,3) - BIG_Px(:,3).*BIG_Py(:,1)));
    h = [sqrt((BIG_Px(:,2) - BIG_Px(:,3)).^2 + (BIG_Py(:,2) - BIG_Py(:,3)).^2)  sqrt((BIG_Px(:,2) - BIG_Px(:,1)).^2 + (BIG_Py(:,2) - BIG_Py(:,1)).^2)  sqrt((BIG_Px(:,3) - BIG_Px(:,1)).^2 + (BIG_Py(:,3) - BIG_Py(:,1)).^2)  max([sqrt((BIG_Px(:,2) - BIG_Px(:,3)).^2 + (BIG_Py(:,2) - BIG_Py(:,3)).^2)  sqrt((BIG_Px(:,2) - BIG_Px(:,1)).^2 + (BIG_Py(:,2) - BIG_Py(:,1)).^2)  sqrt((BIG_Px(:,3) - BIG_Px(:,1)).^2 + (BIG_Py(:,3) - BIG_Py(:,1)).^2)],[],2)];
    Qk = 1./((sqrt(3)/6)*h(:,4))./(A./(0.5*(sum(h(:,1:3),2))));   
    
    [Mx, W] = max(Qk.*[count == 3]'.*in.*Col.*[abs(1-ang/pi) > tolBG]);

    Pi = BIG_PID(W,2); Pip1 = BIG_PID(W,3); Pim1 = BIG_PID(W,1);

    % create a triangle
    TRI.vertices(end+1,:) = [POLYGON.pos(Pim1) POLYGON.pos(Pi) POLYGON.pos(Pip1)];
    TRI.pos(end+1) = length(TRI.pos) + 1;
    TRI.ID(end+1) = length(TRI.ID) + offsetE + 1;

    % Append the list of triangles that a vertex touches
    VER.triangle(POLYGON.pos(Pim1),find(VER.triangle(POLYGON.pos(Pim1),:) == 0,1)) = TRI.pos(end);
    VER.triangle(POLYGON.pos(Pi),find(VER.triangle(POLYGON.pos(Pi),:) == 0,1)) = TRI.pos(end);
    VER.triangle(POLYGON.pos(Pip1),find(VER.triangle(POLYGON.pos(Pip1),:) == 0,1)) = TRI.pos(end);

    % SET 1
    % define edges & assign to element list
    EDGE.vertices(end+1,:) = [min([POLYGON.pos(Pi) POLYGON.pos(Pip1)]) max([POLYGON.pos(Pi) POLYGON.pos(Pip1)])]; EDGE.pos(end+1) = length(EDGE.pos)+1;
    EDGE.triangle(size(EDGE.vertices,1),1) = TRI.pos(end);
    TRI.edges(TRI.pos(end),find(TRI.edges(TRI.pos(end),:) == 0,1)) = EDGE.pos(end);

    % check relative position of vertices Pi and Pip1
    if abs(POLYGON.pos(Pi) - POLYGON.pos(Pip1)) == 1 % vertices are consecutive in order - hence on the boundary
        EDGE.boundary(end+1) = 1;
    elseif POLYGON.pos(Pi) == VER.pos(1) && POLYGON.pos(Pip1) == VER.pos(end) % vertices are on the 'link', also beside each other on the boundary
        EDGE.boundary(end+1) = 1;
    elseif POLYGON.pos(Pi) == VER.pos(end) && POLYGON.pos(Pip1) == VER.pos(1) % vertices are on the 'link', also beside each other on the boundary
        EDGE.boundary(end+1) = 1;
    else % vertices are not on the boundary
        EDGE.boundary(end+1) = 0;
    end

    % assign edge to previously-created element
    if EDGE.boundary(end) == 0
        int = intersect(VER.triangle(EDGE.vertices(end,1),:),VER.triangle(EDGE.vertices(end,2),:));
        if length(int) > 1
            TRI.edges(int(2),find(TRI.edges(int(2),:) == 0,1)) = EDGE.pos(end);
        end
        EDGE.triangle(size(EDGE.vertices,1),2) = int(2);
        EDGE.triangle(size(EDGE.vertices,1),:) = sort(EDGE.triangle(size(EDGE.vertices,1),:));
    end

    % SET 2
    % define edges & assign to element list
    EDGE.vertices(end+1,:) = [min([POLYGON.pos(Pi) POLYGON.pos(Pim1)]) max([POLYGON.pos(Pi) POLYGON.pos(Pim1)])]; EDGE.pos(end+1) = length(EDGE.pos)+1;
    EDGE.triangle(size(EDGE.vertices,1),1) = TRI.pos(end);
    TRI.edges(TRI.pos(end),find(TRI.edges(TRI.pos(end),:) == 0,1)) = EDGE.pos(end);

    % Check relative position of vertices Pi and Pim1
    if abs(POLYGON.pos(Pi) - POLYGON.pos(Pim1)) == 1 % vertices are consecutive in order - hence on the boundary
        EDGE.boundary(end+1) = 1;
    elseif POLYGON.pos(Pi) == VER.pos(1) && POLYGON.pos(Pim1) == VER.pos(end) % vertices are on the 'link', also beside each other on the boundary
        EDGE.boundary(end+1) = 1;
    elseif POLYGON.pos(Pi) == VER.pos(end) && POLYGON.pos(Pim1) == VER.pos(1) % vertices are on the 'link', also beside each other on the boundary
        EDGE.boundary(end+1) = 1;
    else % vertices are not on the boundary
        EDGE.boundary(end+1) = 0;
    end

    % assign edge to previously-created element
    if EDGE.boundary(end) == 0
        int = intersect(VER.triangle(EDGE.vertices(end,1),:),VER.triangle(EDGE.vertices(end,2),:));
        if length(int) > 1
            TRI.edges(int(2),find(TRI.edges(int(2),:) == 0,1)) = EDGE.pos(end);
        end
        EDGE.triangle(size(EDGE.vertices,1),2) = int(2);
        EDGE.triangle(size(EDGE.vertices,1),:) = sort(EDGE.triangle(size(EDGE.vertices,1),:));
    end

    % updating the boundary polygon
    POLYGON.x(Pi) = [];
    POLYGON.y(Pi) = [];
    POLYGON.pos(Pi) = [];
    POLYGON.ID(Pi) = [];
    
end

if length(POLYGON.ID) == 3
    % define final triangle
    TRI.vertices(end+1,:) = [POLYGON.pos(1) POLYGON.pos(2) POLYGON.pos(3)];
    TRI.pos(end+1) = length(TRI.pos) + 1;
    TRI.ID(end+1) = max(TRI.ID) + 1;

    % SET 1
    % define final edges
    EDGE.vertices(end+1,:) = [POLYGON.pos(1) POLYGON.pos(2)]; EDGE.pos(end+1) = length(EDGE.pos)+1;
    EDGE.triangle(size(EDGE.vertices,1),1) = TRI.pos(end);
    TRI.edges(TRI.pos(end),find(TRI.edges(TRI.pos(end),:) == 0,1)) = EDGE.pos(end);

    % check relative position of vertices 1 and 2
    if abs(POLYGON.pos(1) - POLYGON.pos(2)) == 1 % vertices are consecutive in order - hence on the boundary
        EDGE.boundary(end+1) = 1;
    elseif POLYGON.pos(1) == VER.pos(1) && POLYGON.pos(2) == VER.pos(end) % vertices are on the 'link', also beside each other on the boundary
        EDGE.boundary(end+1) = 1;
    elseif POLYGON.pos(1) == VER.pos(end) && POLYGON.pos(2) == VER.pos(1) % vertices are on the 'link', also beside each other on the boundary
        EDGE.boundary(end+1) = 1;
    else % vertices are not on the boundary
        EDGE.boundary(end+1) = 0;
    end

    % assign edge to previously-created element
    if EDGE.boundary(end) == 0
        int = intersect(VER.triangle(EDGE.vertices(end,1),:),VER.triangle(EDGE.vertices(end,2),:));
        if length(int) > 1
            if int(1) == 0
                TRI.edges(int(2),find(TRI.edges(int(2),:) == 0,1)) = EDGE.pos(end);
                EDGE.triangle(size(EDGE.vertices,1),2) = int(2);
            else
                TRI.edges(int(1),find(TRI.edges(int(1),:) == 0,1)) = EDGE.pos(end);
                EDGE.triangle(size(EDGE.vertices,1),2) = int(1);
            end
        end
        EDGE.triangle(size(EDGE.vertices,1),:) = sort(EDGE.triangle(size(EDGE.vertices,1),:));
    end

    % SET 2
    % define final edges
    EDGE.vertices(end+1,:) = [POLYGON.pos(2) POLYGON.pos(3)]; EDGE.pos(end+1) = length(EDGE.pos)+1;
    EDGE.triangle(size(EDGE.vertices,1),1) = TRI.pos(end);
    TRI.edges(TRI.pos(end),find(TRI.edges(TRI.pos(end),:) == 0,1)) = EDGE.pos(end);

    % check relative position of vertices 2 and 3
    if abs(POLYGON.pos(2) - POLYGON.pos(3)) == 1 % vertices are consecutive in order - hence on the boundary
        EDGE.boundary(end+1) = 1;
    elseif POLYGON.pos(2) == VER.pos(1) && POLYGON.pos(3) == VER.pos(end) % vertices are on the 'link', also beside each other on the boundary
        EDGE.boundary(end+1) = 1;
    elseif POLYGON.pos(2) == VER.pos(end) && POLYGON.pos(3) == VER.pos(1) % vertices are on the 'link', also beside each other on the boundary
        EDGE.boundary(end+1) = 1;
    else % vertices are not on the boundary
        EDGE.boundary(end+1) = 0;
    end
    
    % assign edge to previously-created element
    if EDGE.boundary(end) == 0
        int = intersect(VER.triangle(EDGE.vertices(end,1),:),VER.triangle(EDGE.vertices(end,2),:));
        if length(int) > 1
            if int(1) == 0
                TRI.edges(int(2),find(TRI.edges(int(2),:) == 0,1)) = EDGE.pos(end);
                EDGE.triangle(size(EDGE.vertices,1),2) = int(2);
            else
                TRI.edges(int(1),find(TRI.edges(int(1),:) == 0,1)) = EDGE.pos(end);
                EDGE.triangle(size(EDGE.vertices,1),2) = int(1);
            end
        end
        EDGE.triangle(size(EDGE.vertices,1),:) = sort(EDGE.triangle(size(EDGE.vertices,1),:));
    end

    % SET 3
    % define final edges
    EDGE.vertices(end+1,:) = [POLYGON.pos(1) POLYGON.pos(3)]; EDGE.pos(end+1) = length(EDGE.pos)+1;
    EDGE.triangle(size(EDGE.vertices,1),1) = TRI.pos(end);
    TRI.edges(TRI.pos(end),find(TRI.edges(TRI.pos(end),:) == 0,1)) = EDGE.pos(end);

    % check relative position of vertices 1 and 3
    if abs(POLYGON.pos(1) - POLYGON.pos(3)) == 1 % vertices are consecutive in order - hence on the boundary
        EDGE.boundary(end+1) = 1;
    elseif POLYGON.pos(1) == VER.pos(1) && POLYGON.pos(3) == VER.pos(end) % vertices are on the 'link', also beside each other on the boundary
        EDGE.boundary(end+1) = 1;
    elseif POLYGON.pos(1) == VER.pos(end) && POLYGON.pos(3) == VER.pos(1) % vertices are on the 'link', also beside each other on the boundary
        EDGE.boundary(end+1) = 1;
    else % vertices are not on the boundary
        EDGE.boundary(end+1) = 0;
    end

    % assign edge to previously-created element
    if EDGE.boundary(end) == 0
        int = intersect(VER.triangle(EDGE.vertices(end,1),:),VER.triangle(EDGE.vertices(end,2),:));
        if length(int) > 1
            if int(1) == 0
                TRI.edges(int(2),find(TRI.edges(int(2),:) == 0,1)) = EDGE.pos(end);
                EDGE.triangle(size(EDGE.vertices,1),2) = int(2);
            else
                TRI.edges(int(1),find(TRI.edges(int(1),:) == 0,1)) = EDGE.pos(end);
                EDGE.triangle(size(EDGE.vertices,1),2) = int(1);
            end
        end
        EDGE.triangle(size(EDGE.vertices,1),:) = sort(EDGE.triangle(size(EDGE.vertices,1),:));
    end
end
VER = rmfield(VER,'triangle');


%%%%%%%%%%%%%%%%% Algorithm 1.2 - Point insertion - page 13 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if order == 1; disp('Algorithm 2/4 - Interior point insertion'); end
if order == 2; disp('Algorithm 2/5 - Interior point insertion'); end

% appending interior vertices
VER.x(end+1:end+length(INTERIOR.x)) = INTERIOR.x;
VER.y(end+1:end+length(INTERIOR.y)) = INTERIOR.y;
VER.boundary(end+1:end+length(INTERIOR.x)) = 0;
VER.pos = [1:length(VER.pos) + length(INTERIOR.x)];
VER.ID(end+1:end+length(INTERIOR.ID)) = INTERIOR.ID;

for I = 1:length(INTERIOR.x); P = I + length(BOUNDARY.x); ver = 0; 
    for E = 1:length(TRI.pos)
        v1 = TRI.vertices(E,1); e1v1 = EDGE.vertices(TRI.edges(E,1),1); e1v2 = EDGE.vertices(TRI.edges(E,1),2);
        v2 = TRI.vertices(E,2); e2v1 = EDGE.vertices(TRI.edges(E,2),1); e2v2 = EDGE.vertices(TRI.edges(E,2),2);
        v3 = TRI.vertices(E,3); e3v1 = EDGE.vertices(TRI.edges(E,3),1); e3v2 = EDGE.vertices(TRI.edges(E,3),2);
        if (e1v1 == v1 && e1v2 == v2) || (e1v1 == v2 && e1v2 == v1); edge1 = TRI.edges(E,1); end
        if (e1v1 == v2 && e1v2 == v3) || (e1v1 == v3 && e1v2 == v2); edge2 = TRI.edges(E,1); end
        if (e1v1 == v3 && e1v2 == v1) || (e1v1 == v1 && e1v2 == v3); edge3 = TRI.edges(E,1); end
        if (e2v1 == v1 && e2v2 == v2) || (e2v1 == v2 && e2v2 == v1); edge1 = TRI.edges(E,2); end
        if (e2v1 == v2 && e2v2 == v3) || (e2v1 == v3 && e2v2 == v2); edge2 = TRI.edges(E,2); end
        if (e2v1 == v3 && e2v2 == v1) || (e2v1 == v1 && e2v2 == v3); edge3 = TRI.edges(E,2); end
        if (e3v1 == v1 && e3v2 == v2) || (e3v1 == v2 && e3v2 == v1); edge1 = TRI.edges(E,3); end
        if (e3v1 == v2 && e3v2 == v3) || (e3v1 == v3 && e3v2 == v2); edge2 = TRI.edges(E,3); end
        if (e3v1 == v3 && e3v2 == v1) || (e3v1 == v1 && e3v2 == v3); edge3 = TRI.edges(E,3); end

        x1 = VER.x(v1); y1 = VER.y(v1); 
        x2 = VER.x(v2); y2 = VER.y(v2); 
        x3 = VER.x(v3); y3 = VER.y(v3);
        xi = INTERIOR.x(I); yi = INTERIOR.y(I);
      
        % checking if point is inside the triangle
        in = EXTERNAL_inpoly([xi yi],[[x1 x2 x3]' [y1 y2 y3]']);
        
        % checking if point is on the edges of the triangle
        Mat12 = [(x1 - x2) (y1 - y2); (xi - x2) (yi - y2)]; Col12 = abs(det(Mat12));
        Mat13 = [(x1 - x3) (y1 - y3); (xi - x3) (yi - y3)]; Col13 = abs(det(Mat13));
        Mat23 = [(x2 - x3) (y2 - y3); (xi - x3) (yi - y3)]; Col23 = abs(det(Mat23));
        if (Col12 < tolNN) || (Col13 < tolNN) || (Col23 < tolNN); on = 1; else on = 0; end        

        % checking if interior point corresponds to an already existing vertex
        if sum((INTERIOR.x(I) == [x1 x2 x3]).*(INTERIOR.y(I) == [y1 y2 y3])) > 0; ver = 1; end

        if in == 1 && on == 1 && ver == 0 % point is on an edge but is not a vertext
            % degenerate edge triangles into two sub-triangles per edge triangle
            vedge = VER.pos(P); xedge = VER.x(P); yedge = VER.y(P);

            % which edge? find out by checking for collinearity
            mat12 = [(x1 - x2) (y1 - y2); (xedge - x2) (yedge - y2)]; col12 = abs(det(mat12));
            mat23 = [(x2 - x3) (y2 - y3); (xedge - x3) (yedge - y3)]; col23 = abs(det(mat23));
            mat13 = [(x1 - x3) (y1 - y3); (xedge - x3) (yedge - y3)]; col13 = abs(det(mat13));
            if col12 < tolNN; edge = edge1; vopp = v3; vother1 = v1; vother2 = v2; eother1 = edge3; eother2 = edge2; end
            if col23 < tolNN; edge = edge2; vopp = v1; vother1 = v2; vother2 = v3; eother1 = edge1; eother2 = edge3; end
            if col13 < tolNN; edge = edge3; vopp = v2; vother1 = v3; vother2 = v1; eother1 = edge2; eother2 = edge1; end

            % second triangle
            wh = EDGE.triangle(edge,find(EDGE.triangle(edge,:) ~= TRI.pos(E)));
            V1 = TRI.vertices(wh,1); E1v1 = EDGE.vertices(TRI.edges(wh,1),1); E1v2 = EDGE.vertices(TRI.edges(wh,1),2);
            V2 = TRI.vertices(wh,2); E2v1 = EDGE.vertices(TRI.edges(wh,2),1); E2v2 = EDGE.vertices(TRI.edges(wh,2),2);
            V3 = TRI.vertices(wh,3); E3v1 = EDGE.vertices(TRI.edges(wh,3),1); E3v2 = EDGE.vertices(TRI.edges(wh,3),2);
            X1 = VER.x(V1); Y1 = VER.y(V1);
            X2 = VER.x(V2); Y2 = VER.y(V2);
            X3 = VER.x(V3); Y3 = VER.y(V3);
            if (E1v1 == V1 && E1v2 == V2) || (E1v1 == V2 && E1v2 == V1); Edge1 = TRI.edges(wh,1); end
            if (E1v1 == V2 && E1v2 == V3) || (E1v1 == V3 && E1v2 == V2); Edge2 = TRI.edges(wh,1); end
            if (E1v1 == V3 && E1v2 == V1) || (E1v1 == V1 && E1v2 == V3); Edge3 = TRI.edges(wh,1); end
            if (E2v1 == V1 && E2v2 == V2) || (E2v1 == V2 && E2v2 == V1); Edge1 = TRI.edges(wh,2); end
            if (E2v1 == V2 && E2v2 == V3) || (E2v1 == V3 && E2v2 == V2); Edge2 = TRI.edges(wh,2); end
            if (E2v1 == V3 && E2v2 == V1) || (E2v1 == V1 && E2v2 == V3); Edge3 = TRI.edges(wh,2); end
            if (E3v1 == V1 && E3v2 == V2) || (E3v1 == V2 && E3v2 == V1); Edge1 = TRI.edges(wh,3); end
            if (E3v1 == V2 && E3v2 == V3) || (E3v1 == V3 && E3v2 == V2); Edge2 = TRI.edges(wh,3); end
            if (E3v1 == V3 && E3v2 == V1) || (E3v1 == V1 && E3v2 == V3); Edge3 = TRI.edges(wh,3); end

            Mat12 = [(X1 - X2) (Y1 - Y2); (xedge - X2) (yedge - Y2)]; Col12 = abs(det(Mat12));
            Mat23 = [(X2 - X3) (Y2 - Y3); (xedge - X3) (yedge - Y3)]; Col23 = abs(det(Mat23));
            Mat13 = [(X1 - X3) (Y1 - Y3); (xedge - X3) (yedge - Y3)]; Col13 = abs(det(Mat13));
            if Col12 < tolNN; Edge = Edge1; Vopp = V3; Vother1 = V1; Vother2 = V2; Eother1 = Edge3; Eother2 = Edge2; end
            if Col23 < tolNN; Edge = Edge2; Vopp = V1; Vother1 = V2; Vother2 = V3; Eother1 = Edge1; Eother2 = Edge3; end
            if Col13 < tolNN; Edge = Edge3; Vopp = V2; Vother1 = V3; Vother2 = V1; Eother1 = Edge2; Eother2 = Edge1; end

            % the first triangle degeneration
            % 1st sub-triangle
            TRI.vertices(end+1,:) = [vother1 vedge vopp];
            TRI.pos(end+1) = length(TRI.pos) + 1; % t1
            TRI.ID(end+1) = max(TRI.ID) + 1;

            % 2nd sub-triangle
            TRI.vertices(end+1,:) = [vopp vedge vother2];
            TRI.pos(end+1) = length(TRI.pos) + 1; % t2
            TRI.ID(end+1) = max(TRI.ID) + 1;

            % create new edges
            EDGE.vertices(end+1,:) = [min([vedge vother1]) max([vedge vother1])]; EDGE.pos(end+1) = EDGE.pos(end) + 1; EDGE.boundary(end+1) = 0; n1 = EDGE.pos(end);
            EDGE.vertices(end+1,:) = [min([vedge vother2]) max([vedge vother2])]; EDGE.pos(end+1) = EDGE.pos(end) + 1; EDGE.boundary(end+1) = 0; n2 = EDGE.pos(end);
            EDGE.vertices(end+1,:) = [min([vedge vopp]) max([vedge vopp])]; EDGE.pos(end+1) = EDGE.pos(end) + 1; EDGE.boundary(end+1) = 0; n3 = EDGE.pos(end);

            % update triangle edges
            TRI.edges(end+1,:) = [n1 n3 eother1];
            TRI.edges(end+1,:) = [n3 n2 eother2];

            % the second triangle denegeration
            % 1st sub-triangle
            TRI.vertices(end+1,:) = [Vother1 vedge Vopp];
            TRI.pos(end+1) = length(TRI.pos) + 1; % t1
            TRI.ID(end+1) = max(TRI.ID) + 1;

            % 2nd sub-triangle
            TRI.vertices(end+1,:) = [Vopp vedge Vother2];
            TRI.pos(end+1) = length(TRI.pos) + 1; % t2
            TRI.ID(end+1) = max(TRI.ID) + 1;

            % create new edge
            EDGE.vertices(end+1,:) = [min([vedge Vopp]) max([vedge Vopp])]; EDGE.pos(end+1) = EDGE.pos(end) + 1; EDGE.boundary(end+1) = 0; N3 = EDGE.pos(end);

            % update triangle edges
            TRI.edges(end+1,:) = [n2 N3 Eother1];
            TRI.edges(end+1,:) = [N3 n1 Eother2];

            % update triangle edge associations
            EDGE.triangle(eother1,find(EDGE.triangle(eother1,:) == TRI.pos(E))) = TRI.pos(end-3);
            EDGE.triangle(eother2,find(EDGE.triangle(eother2,:) == TRI.pos(E))) = TRI.pos(end-2);
            EDGE.triangle(Eother1,find(EDGE.triangle(Eother1,:) == TRI.pos(wh))) = TRI.pos(end-1);
            EDGE.triangle(Eother2,find(EDGE.triangle(Eother2,:) == TRI.pos(wh))) = TRI.pos(end);

            EDGE.triangle(end+1,:) = [min([TRI.pos(end-3) TRI.pos(end)]) max([TRI.pos(end-3) TRI.pos(end)])];
            EDGE.triangle(end+1,:) = [min([TRI.pos(end-2) TRI.pos(end-1)]) max([TRI.pos(end-2) TRI.pos(end-1)])];
            EDGE.triangle(end+1,:) = [min([TRI.pos(end-3) TRI.pos(end-2)]) max([TRI.pos(end-3) TRI.pos(end-2)])];
            EDGE.triangle(end+1,:) = [min([TRI.pos(end-1) TRI.pos(end)]) max([TRI.pos(end-1) TRI.pos(end)])];

            % delete the original triangle and edge
            TRI.vertices([TRI.pos(E) wh],:) = [];
            TRI.edges([TRI.pos(E) wh],:) = [];
            TRI.pos([TRI.pos(E) wh]) = [];
            TRI.ID([TRI.pos(E) wh]) = [];

            EDGE.vertices(edge,:) = [];
            EDGE.triangle(edge,:) = [];
            EDGE.boundary(edge) = [];
            EDGE.pos(edge) = [];

            for p = 1:length(TRI.pos)
                old = TRI.pos(p); TRI.pos(p) = p;
                EDGE.triangle(find(EDGE.triangle == old)) = p;
            end

            for e = 1:length(EDGE.pos)
                old = EDGE.pos(e); EDGE.pos(e) = e;
                TRI.edges(find(TRI.edges == old)) = e;
            end
            
        elseif in == 1 && on == 0 && ver == 0 % point is fully inside the triangle and is not a vertex
            % degenerate triangle into three sub-triangles
            vint = VER.pos(P);

            % 1st sub-triangle
            TRI.vertices(end+1,:) = [v1 v2 vint];
            TRI.pos(end+1) = length(TRI.pos) + 1;
            TRI.ID(end+1) = max(TRI.ID) + 1;

            % 2nd sub-triangle
            TRI.vertices(end+1,:) = [v2 v3 vint];
            TRI.pos(end+1) = length(TRI.pos) + 1;
            TRI.ID(end+1) = max(TRI.ID) + 1;

            % 3rd sub-triangle
            TRI.vertices(end+1,:) = [v3 v1 vint];
            TRI.pos(end+1) = length(TRI.pos) + 1;
            TRI.ID(end+1) = max(TRI.ID) +  1;

            % adding extra interior edges and updating edge adjacency data
            EDGE.vertices(end+1,:) = [v2 vint]; EDGE.pos(end+1) = EDGE.pos(end) + 1; EDGE.boundary(end+1) = 0;
            EDGE.vertices(end+1,:) = [v3 vint]; EDGE.pos(end+1) = EDGE.pos(end) + 1; EDGE.boundary(end+1) = 0;
            EDGE.vertices(end+1,:) = [v1 vint]; EDGE.pos(end+1) = EDGE.pos(end) + 1; EDGE.boundary(end+1) = 0;

            EDGE.triangle(end+1,:) = [TRI.pos(end-2) TRI.pos(end-1)];
            EDGE.triangle(end+1,:) = [TRI.pos(end-1) TRI.pos(end)];
            EDGE.triangle(end+1,:) = [TRI.pos(end)   TRI.pos(end-2)];

            EDGE.triangle(edge1,find(EDGE.triangle(edge1,:) == TRI.pos(E))) = TRI.pos(end-2);
            EDGE.triangle(edge2,find(EDGE.triangle(edge2,:) == TRI.pos(E))) = TRI.pos(end-1);
            EDGE.triangle(edge3,find(EDGE.triangle(edge3,:) == TRI.pos(E))) = TRI.pos(end);

            TRI.edges(end+1,:) = [edge1 EDGE.pos(end-2) EDGE.pos(end)];
            TRI.edges(end+1,:) = [edge2 EDGE.pos(end-1) EDGE.pos(end-2)];
            TRI.edges(end+1,:) = [edge3 EDGE.pos(end)   EDGE.pos(end-1)];

            % delete the original triangle
            TRI.vertices(E,:) = [];
            TRI.edges(E,:) = [];
            TRI.pos(E) = [];
            TRI.ID(E) = [];

            for p = 1:length(TRI.pos)
                old = TRI.pos(p); TRI.pos(p) = p;
                EDGE.triangle(find(EDGE.triangle == old)) = p;
            end
         end
    end
end


%%%%%%%%%%%%%%%%% Algorithms 1.3 & 3.2 - Edge swapping/Circumcentre/LOP - pages 15 & 66 %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if order == 1; disp('Algorithm 3/4 - Edge swapping and local optimisation'); end
if order == 2; disp('Algorithm 3/5 - Edge swapping and local optimisation'); end

for G = 1:length(EDGE.boundary)
    if EDGE.boundary(G) == 0 && EDGE.boundary(G) ~= -1 % interior edge
        % identifying adjacent triangles and relative vertices
        T1 = EDGE.triangle(G,1); v1t1 = TRI.vertices(T1,1); v2t1 = TRI.vertices(T1,2); v3t1 = TRI.vertices(T1,3);
        T2 = EDGE.triangle(G,2); v1t2 = TRI.vertices(T2,1); v2t2 = TRI.vertices(T2,2); v3t2 = TRI.vertices(T2,3);
        a = [v1t1 v2t1 v3t1]; b = [v1t2 v2t2 v3t2]; d = intersect(a,b);
        V2 = a(a == d(1));
        V3 = a(a == d(2));
        a(a == d(1)) = []; a(a == d(2)) = []; b(b == d(1)) = []; b(b == d(2)) = [];
        V1 = a; V4 = b;

        % identifying relative edges
        tri1 = find(TRI.edges(T1,:) ~= EDGE.pos(G)); TRI1 = TRI.edges(T1,tri1);
        tri2 = find(TRI.edges(T2,:) ~= EDGE.pos(G)); TRI2 = TRI.edges(T2,tri2);
        e1t1 = TRI1(1); e2t1 = TRI1(2); e1t2 = TRI2(1); e2t2 = TRI2(2);
        if sum(EDGE.vertices(e1t1,:) == [V1 V2]) == 2 || sum(EDGE.vertices(e1t1,:) == [V2 V1]) == 2; T1E1 = e1t1; T1E2 = e2t1; end
        if sum(EDGE.vertices(e1t1,:) == [V1 V3]) == 2 || sum(EDGE.vertices(e1t1,:) == [V3 V1]) == 2; T1E1 = e2t1; T1E2 = e1t1; end
        if sum(EDGE.vertices(e1t2,:) == [V2 V4]) == 2 || sum(EDGE.vertices(e1t2,:) == [V4 V2]) == 2; T2E1 = e1t2; T2E2 = e2t2; end
        if sum(EDGE.vertices(e1t2,:) == [V3 V4]) == 2 || sum(EDGE.vertices(e1t2,:) == [V4 V3]) == 2; T2E1 = e2t2; T2E2 = e1t2; end

        % nodal coordinates and the circumcentre test (checking if convex quadrilateral)
        x1 = VER.x(V1); x2 = VER.x(V2); x3 = VER.x(V3); x4 = VER.x(V4);
        y1 = VER.y(V1); y2 = VER.y(V2); y3 = VER.y(V3); y4 = VER.y(V4);
        a213 = [x1 - x2, y1 - y2]; a243 = [x4 - x2, y4 - y2];
        b213 = [x1 - x3, y1 - y3]; b243 = [x4 - x3, y4 - y3];
        ang213 = acos(dot(a213,b213)/(norm(a213)*norm(b213)));
        ang243 = acos(dot(a243,b243)/(norm(a243)*norm(b243)));
        if ang213 + ang243 > pi; D = 1; else D = 0; end

        % convex hull number of points verification
        %K = convhull([x1 x2 x3 x4],[y1 y2 y3 y4],{'Qt','Pp'}); K(end) = [];
        K = convhull([x1 x2 x3 x4],[y1 y2 y3 y4]); K(end) = [];
        if length(K) == 4; C = 1; else C = 0; end

        % check for collinearity of vertices - such edges should not be swapped
        Mat124 = [(x1 - x4) (y1 - y4); (x2 - x4) (y2 - y4)]; Col124 = abs(det(Mat124));
        Mat134 = [(x1 - x4) (y1 - y4); (x3 - x4) (y3 - y4)]; Col134 = abs(det(Mat134));

        if (D == 1) && (Col124 > tolBG) && (Col134 > tolBG) && (C == 1);
            % updating triangle vertices and edge associations
            TRI.vertices(T1,:) = [V1 V2 V4];
            TRI.vertices(T2,:) = [V1 V3 V4];
            TRI.edges(T1,:) = [T1E1 T2E1 EDGE.pos(G)];
            TRI.edges(T2,:) = [T1E2 T2E2 EDGE.pos(G)];
            EDGE.vertices(G,:) = [V1 V4];
            EDGE.boundary(G) = -1;
            EDGE.triangle(T1E2,find(EDGE.triangle(T1E2,:) == T1)) = T2;
            EDGE.triangle(T2E1,find(EDGE.triangle(T2E1,:) == T2)) = T1;
        else % edge is already locally optimal - mark it but do nothing more with it
            EDGE.boundary(G) = -1;
        end
    end
end


%%%%%%%%%%%%%%%%% Addition of mid-edge nodes for second order triangulation %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if order == 2; first = 1;
    disp('Algorithm 4/5 - Converting to second order triangulation');
    TRI.vertices(:,4:6) = zeros(size(TRI.vertices));
    EDGE.vertices(:,3) = zeros(size(EDGE.vertices(:,1)));
    for E = 1:length(EDGE.pos)
        v1 = EDGE.vertices(E,1); x1 = VER.x(v1); y1 = VER.y(v1);
        v2 = EDGE.vertices(E,2); x2 = VER.x(v2); y2 = VER.y(v2);       
        
        % creating extre mid-edge vertex
        VER.x(end+1) = 0.5*(x1 + x2);
        VER.y(end+1) = 0.5*(y1 + y2);
        if (VER.boundary(EDGE.vertices(E,1)) == 1) && (VER.boundary(EDGE.vertices(E,2)) == 1) 
            VER.boundary(end+1) = 1; 
        else
            VER.boundary(end+1) = 0;
        end
        VER.pos(end+1) = max(VER.pos) + 1;
        if first == 1; VER.ID(end+1) = offsetN; first = 0; else VER.ID(end+1) = max(VER.ID) + 1; end
        EDGE.vertices(E,3) = VER.pos(end);
        
        % appending triangle vertices
        T1 = EDGE.triangle(E,1);
        T2 = EDGE.triangle(E,2);
        if T1 > 0
            a = find(TRI.vertices(T1,:) == v1);
            b = find(TRI.vertices(T1,:) == v2);
            if ((a == 1) && (b == 2)) || ((a == 2) && (b == 1)); TRI.vertices(T1,4) = VER.pos(end); end
            if ((a == 2) && (b == 3)) || ((a == 3) && (b == 2)); TRI.vertices(T1,5) = VER.pos(end); end
            if ((a == 1) && (b == 3)) || ((a == 3) && (b == 1)); TRI.vertices(T1,6) = VER.pos(end); end          
        end
        if T2 > 0
            a = find(TRI.vertices(T2,:) == v1);
            b = find(TRI.vertices(T2,:) == v2);
            if ((a == 1) && (b == 2)) || ((a == 2) && (b == 1)); TRI.vertices(T2,4) = VER.pos(end); end
            if ((a == 2) && (b == 3)) || ((a == 3) && (b == 2)); TRI.vertices(T2,5) = VER.pos(end); end
            if ((a == 1) && (b == 3)) || ((a == 3) && (b == 1)); TRI.vertices(T2,6) = VER.pos(end); end          
        end        
        no = find(TRI.edges);
    end
end


%%%%%%%%%%%%%%%%% Boundary element edge normalisation %%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if order == 1; disp('Algorithm 4/4 - Boundary element edge consistency check'); end
if order == 2; disp('Algorithm 5/5 - Boundary element edge consistency check'); end

boundaries = find(EDGE.boundary == 1); TRIold = TRI.vertices;
TRI.vertices = VER.ID(TRI.vertices); TRI_boundary = [];
EDGE.vertices = VER.ID(EDGE.vertices); VER_boundary = [];
for E = 1:length(boundaries)
    edge = boundaries(E); T = EDGE.triangle(edge,1); if T == 0; T = EDGE.triangle(edge,2); end
    v1 = EDGE.vertices(edge,1); a = find(TRI.vertices(T,:) == v1);
    v2 = EDGE.vertices(edge,2); b = find(TRI.vertices(T,:) == v2);

    if ((a == 1) && (b == 2)) || ((a == 2) && (b == 1)) % case 1 rearrangement
        V1 = TRI.vertices(T,1); V2 = TRI.vertices(T,2); V3 = TRI.vertices(T,3);
        if order == 2
            V4 = TRI.vertices(T,4); V5 = TRI.vertices(T,5); V6 = TRI.vertices(T,6);
            TRI.vertices(T,:) = [V2 V3 V1 V5 V6 V4];
            VER_boundary(end+1) = V2; VER_boundary(end+1) = V1; VER_boundary(end+1) = V4;
        else
            TRI.vertices(T,:) = [V2 V3 V1];
            VER_boundary(end+1) = V2; VER_boundary(end+1) = V1; 
        end
    elseif ((a == 2) && (b == 3)) || ((a == 3) && (b == 2)); % case 2 rearrangement
        V1 = TRI.vertices(T,1); V2 = TRI.vertices(T,2); V3 = TRI.vertices(T,3);
        if order == 2
            V4 = TRI.vertices(T,4); V5 = TRI.vertices(T,5); V6 = TRI.vertices(T,6);
            TRI.vertices(T,:) = [V3 V1 V2 V6 V4 V5];
            VER_boundary(end+1) = V3; VER_boundary(end+1) = V2; VER_boundary(end+1) = V5;
        else
            TRI.vertices(T,:) = [V3 V1 V2];
            VER_boundary(end+1) = V3; VER_boundary(end+1) = V2; 
        end
    else % case 3 is neutral - no rearrangement
        VER_boundary(end+1) = TRI.vertices(T,1); VER_boundary(end+1) = TRI.vertices(T,3);
        if order == 2; VER_boundary(end+1) = TRI.vertices(T,6); end
    end
    TRI_boundary(end+1) = T;
end
clear('EDGE');
VER = rmfield(VER,'boundary'); VER = rmfield(VER,'pos');
TRI = rmfield(TRI,'edges'); TRI = rmfield(TRI,'pos'); TRI_boundary = TRI.ID(TRI_boundary);
VER_boundary = unique(VER_boundary); TRI_boundary = unique(TRI_boundary);
VER.x = VER.x/factX; VER.y = VER.y/factY;