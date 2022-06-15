function [NODE,ELEMENT_RECT,OP] = AJSABQ_helical_phoenix(PP,RSL,OP,NODE,ELEMENT_RECT,GEO,MATERIAL)
% 1st part of the MPI Parallelised helical meshing algorithm
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:24 (previously 01/08/12 - 19:52)

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

disp(['Main meshing (helical shell elements): Rank ',num2str(PP.rank)]);
master = PP.master; slave = PP.slave; rank = PP.rank; comm = PP.comm; CPUs = PP.CPUs; conformal = GEO.helix.conformal;
Pitch = GEO.helix.pitch; np = GEO.helix.no_p; alpha = GEO.helix.incline*pi/180; Alpha = GEO.helix.incline;
rho = 2*pi*OP.R*sin(alpha);

% Defining functions to express both the left and right diagonal lines
fL_y = @(x,PP,m)   PP - x*m;
fL_x =  @(y,PP,m)  (PP - y)/m;
fR_y = @(x,xmax,PP,m)   m*(xmax - x);
fR_x =  @(y,xmax,PP,m)  xmax - y/m;
arc = @(phi,P,R) phi*sqrt(R^2 + (P/(2*pi))^2);

if conformal == 0
    PPP = Pitch; % pitch - zeta, phi plane
    x_max = 2*(np+1)*pi; % max phi
    m = Pitch/(2*pi); % Local slope angle in zeta - phi plane (beta)
elseif conformal == 1
    PPP = 2*pi*OP.R*sin(alpha); % rho - eta, sigma plane
    x_max = arc(2*(np)*pi,Pitch,OP.R) + 2*pi*OP.R*cos(alpha); % max sigma
    m = tan(alpha); % Local slope angle in eta - sigma plane (alpha)
end

tolS = OP.Tols.tolS; % Side tolerance
tolH = OP.Tols.tolH; % Vertical tolerance
xoff = OP.Tols.sideE*abs(RSL.t2(1) - RSL.t1(1))/RSL.T(1); OP.xoff = xoff; % Horizontal offset
P = struct('Lefts',struct('x',[],'y',[],'ID',[]),'Lefts2',struct('x',[],'y',[],'ID',[]),...
    'Rights',struct('x',[],'y',[],'ID',[]),'Rights2',struct('x',[],'y',[],'ID',[]),...
    'Tops',struct('x',[],'y',[],'ID',[]),'Tops2',struct('x',[],'y',[],'ID',[]),...
    'Ints',struct('x',[],'y',[],'ID',[]),'Ints2',struct('x',[],'y',[],'ID',[]),...
    'Mids',struct('x',[],'y',[],'ID',[]),'Mids2',struct('x',[],'y',[],'ID',[]),...
    'Bots',struct('x',[],'y',[],'ID',[]),'Bots2',struct('x',[],'y',[],'ID',[]),...
    'LTop',struct('x',[],'y',[],'ID',[]),'LBot',struct('x',[],'y',[],'ID',[]),...
    'RTop',struct('x',[],'y',[],'ID',[]),'RBot',struct('x',[],'y',[],'ID',[]),...
    'Ls',struct('x',[],'y',[]),'Rs',struct('x',[],'y',[])); % Node sets for the triangulation
DEL = struct('N',[],'E',[]); % Marking of nodes/elements to be deleted later

i1 = OP.i1; i2 = OP.i2; e1 = OP.e1; e2 = OP.e2; div = OP.div; Ymin = OP.Hmin; Ymax = OP.Hmax; Ttot = max(RSL.t2); Gmax = length(RSL.Region);
P.LTop.x = fL_x(PPP,PPP,m); P.LBot.x = fL_x(0,PPP,m); P.RTop.x = fR_x(PPP,x_max,PPP,m); P.RBot.x = fR_x(0,x_max,PPP,m);
P.LTop.y = PPP; P.LBot.y = 0; P.RTop.y = PPP; P.RBot.y = 0;

if isempty(find(master == rank,1)) == 0; N_top = []; N_bot = []; end

for G = 1:length(RSL.Region)
    for D = 1:CPUs
        if isempty(find(master == rank,1)) == 0; D = 1; end
        if isempty(find(slave == rank,1)) == 0; D = PP.rank + 1;
            NODE = struct('index',[],'x',[],'y',[],'X',[],'Y',[],'Z',[],'case',[]);
            ELEMENT_RECT = struct('index',[],'nodes',[],'type',[],'case',[]);
        end
    end

    Imin = i1(G,D); Imax = i2(G,D); ymin = Ymin(G,D); ymax = Ymax(G,D); Emin = e1(G,D); Emax = e2(G,D); I = Imin - 1; E = Emin - 1; RTcurrent = RSL.Ttot+1; dB = 0; dT = 0; DX = []; first = 1; First = 1;
    if G == 1 && D == 1; Yno = abs(ymax - ymin)/div(G,D); else Yno = abs(ymax - ymin)/(div(G,D)-1); end

    for Y = ymin:Yno:ymax
        for A = 1:length(RSL.t1)
            xmax = RSL.t2(A); xmin = RSL.t1(A); Xno = abs(xmax - xmin)/RSL.T(A);
            if A == 1; xmin = xmin; else xmin = xmin + Xno; end
            for X = xmin:Xno:xmax; I = I + 1;
                if Y < tolH;
                    if first == 1; xprev = X; first = 0;
                    else dT = dT + 1; DX(dT) = X - xprev;
                    end
                end

                if Y < tolH && X >= fL_x(Y,PPP,m) - tolS
                    if First == 1; Xtemp = X; Xprev = X; First = 0;
                    else
                        dB = dB + 1; Xtemp = Xprev + DX(dB);
                    end
                else
                    Xtemp = X;
                end

                % Spiral Welds (if present)
                if GEO.feature.H_Weld.Toggle == 1 % formulation adapted from Rotter and Teng, 1989
                    t = GEO.helix.thick;
                    nu = MATERIAL.v;
                    d0 = GEO.feature.H_Weld.Ampt;
                    lambda = pi*sqrt(OP.R*t)*(3*(1-nu*nu))^0.25;
                    if GEO.feature.H_Weld.Type == 'A'
                        N_R = inline('R-d0*exp(-pi*zeta/lambda)*(cos(pi*zeta/lambda) + sin(pi*zeta/lambda))','R','d0','zeta','lambda');
                    elseif GEO.feature.H_Weld.Type == 'B'
                        N_R = inline('R-d0*exp(-pi*zeta/lambda)*(cos(pi*zeta/lambda))','R','d0','zeta','lambda');
                    end
                else
                    d0 = 0; lambda = 0; N_R = inline('R','R','d0','zeta','lambda');
                end

                % Radial coordinate
                NR = N_R(OP.R,d0,abs(Y-0.5*PPP),lambda);

                %%% ASSIGNING NODES
                if isempty(find(master == rank,1)) == 0
                    NODE.index(I) = I;

                    % 2D coordinates
                    NODE.x(I) = Xtemp; NODE.y(I) = Y; NODE.case(I) = 2;
                    if X < fL_x(Y,PPP,m) + tolS + xoff; NODE.case(I) = 1; DEL.N(end+1) = I; end
                    if X > fR_x(Y,x_max,PPP,m) - tolS - xoff; NODE.case(I) = 3; DEL.N(end+1) = I; end
                    if (X <= fL_x(Y,PPP,m) + tolS + xoff) && (X >= fL_x(Y,PPP,m) - tolS + xoff); NODE.case(I) = 1; DEL.N(end+1) = I; end
                    if (X <= fR_x(Y,x_max,PPP,m) - tolS - xoff) && (X >= fR_x(Y,x_max,PPP,m) + tolS - xoff); NODE.case(I) = 3; DEL.N(end+1) = I; end                    
                    if (X <= fL_x(Y,PPP,m) + (2/3)*xoff - tolS) && (X >= fL_x(Y,PPP,m) + (1/3)*xoff + tolS) && (Y >= 0.0 + tolS && Y <= PPP - tolS); P.Ls.x(end+1) = X; P.Ls.y(end+1) = Y; end
                    if (X <= fR_x(Y,x_max,PPP,m) - (1/3)*xoff - tolS) && (X >= fR_x(Y,x_max,PPP,m) - (2/3)*xoff + tolS) && (Y >= 0.0 + tolS && Y <= PPP - tolS); P.Rs.x(end+1) = X; P.Rs.y(end+1) = Y; end
                    
                    % 3D coordinates
                    [NODE.X(I), NODE.Y(I), NODE.Z(I)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODE.x(I),NODE.y(I),NR);
                end
                if isempty(find(slave == rank,1)) == 0
                    NODE.index(I+1-Imin) = I;

                    % 2D coordinates
                    NODE.x(I+1-Imin) = Xtemp; NODE.y(I+1-Imin) = Y; NODE.case(I+1-Imin) = 2;
                    if X < fL_x(Y,PPP,m) + tolS + xoff; NODE.case(I+1-Imin) = 1; DEL.N(end+1) = I; end
                    if X > fR_x(Y,x_max,PPP,m) - tolS - xoff; NODE.case(I+1-Imin) = 3; DEL.N(end+1) = I; end
                    if (X <= fL_x(Y,PPP,m) + tolS + xoff) && (X >= fL_x(Y,PPP,m) - tolS + xoff); NODE.case(I+1-Imin) = 1; DEL.N(end+1) = I; end
                    if (X <= fR_x(Y,x_max,PPP,m) - tolS - xoff) && (X >= fR_x(Y,x_max,PPP,m) + tolS - xoff); NODE.case(I+1-Imin) = 3; DEL.N(end+1) = I; end                    
                    if (X <= fL_x(Y,PPP,m) + (2/3)*xoff - tolS) && (X >= fL_x(Y,PPP,m) + (1/3)*xoff + tolS) && (Y >= 0.0 + tolS && Y <= PPP - tolS); P.Ls.x(end+1) = X; P.Ls.y(end+1) = Y; end
                    if (X <= fR_x(Y,x_max,PPP,m) - (1/3)*xoff - tolS) && (X >= fR_x(Y,x_max,PPP,m) - (2/3)*xoff + tolS) && (Y >= 0.0 + tolS && Y <= PPP - tolS); P.Rs.x(end+1) = X; P.Rs.y(end+1) = Y; end
                    
                    % 3D coordinates
                    [NODE.X(I+1-Imin), NODE.Y(I+1-Imin), NODE.Z(I+1-Imin)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODE.x(I+1-Imin),NODE.y(I+1-Imin),NR);
                end
            end
        end
    end

    %%% CALIBRATION OF NODES
    if isempty(find(slave == rank,1)) == 0; MPI_Send(master,4+G+rank,comm,NODE); end
    if isempty(find(master == rank,1)) == 0 && CPUs > 1
        for S = 1:length(slave)
            [SlaveNodes] = MPI_Recv(slave(S),4+G+slave(S),comm);
            NODE.index(SlaveNodes.index) = SlaveNodes.index;
            NODE.x(SlaveNodes.index) = SlaveNodes.x;
            NODE.y(SlaveNodes.index) = SlaveNodes.y;
            NODE.X(SlaveNodes.index) = SlaveNodes.X;
            NODE.Y(SlaveNodes.index) = SlaveNodes.Y;
            NODE.Z(SlaveNodes.index) = SlaveNodes.Z;
            NODE.case(SlaveNodes.index) = SlaveNodes.case;
        end

        % Broadcasting updated node struct (so far) to other threads
        MPI_Bcast(master,1000+G,comm,NODE);
    end
    if isempty(find(slave == rank,1)) == 0 && CPUs > 1; [NODE] = MPI_Recv(master,1000+G,comm); end

    if G == 1 && CPUs == 1; BaseZ = 1; elseif G > 1 && CPUs == 1; BaseZ = sum(RSL.Z(1:G-1)) + 1; end
    if CPUs == 1; TopZ1 = BaseZ + RSL.Z(G) - 1; TopZ2 = BaseZ + RSL.Z(G) - 2; end
    if G == 1 && CPUs > 1; BaseZ = 1 + sum(div(G,1:rank)); elseif G > 1 && CPUs > 1; BaseZ = sum(RSL.Z(1:G-1)) + sum(div(G,1:rank)) + 1; end
    if CPUs > 1; TopZ1 = BaseZ + div(G,rank+1) - 1; TopZ2 = BaseZ + div(G,rank+1) - 2; end
    if isempty(find(master == rank,1)) == 0; J = E; end
    if isempty(find(slave == rank,1)) == 0; J = E + 1 - Emin; end

    %%% ELEMENT GENERATION
    if GEO.helix.element == 'a' || GEO.helix.element == 'b' || GEO.helix.element == 'c'; % If it is a 4-node element
        for Zn = BaseZ:TopZ1 % For every row but one of vertical nodes
            for Tn = 1:RSL.Ttot; E = E + 1; J = J + 1; % For every column but one of circumferential nodes
                if Tn == 1; first = 1; end
                ELEMENT_RECT.index(J) = E; ELEMENT_RECT.type(J) = GEO.helix.element; ELEMENT_RECT.case(J) = 0;
                ELEMENT_RECT.nodes(J,1) = Tn + (Zn-1)*RTcurrent;
                ELEMENT_RECT.nodes(J,2) = Tn + 1 + (Zn-1)*RTcurrent;
                ELEMENT_RECT.nodes(J,3) = Tn + 1 + Zn*RTcurrent;
                ELEMENT_RECT.nodes(J,4) = Tn + Zn*RTcurrent;

                if Tn == 1
                    dx = abs(NODE.x(ELEMENT_RECT.nodes(J,1)) - NODE.x(ELEMENT_RECT.nodes(J,2)));
                    dy = abs(NODE.y(ELEMENT_RECT.nodes(J,1)) - NODE.y(ELEMENT_RECT.nodes(J,4)));
                    [X1 Y1 Z1] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODE.x(ELEMENT_RECT.nodes(J,1)),NODE.y(ELEMENT_RECT.nodes(J,1)),N_R(OP.R,d0,abs(NODE.y(ELEMENT_RECT.nodes(J,1))-0.5*PPP),lambda));
                    [X2 Y2 Z2] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODE.x(ELEMENT_RECT.nodes(J,2)),NODE.y(ELEMENT_RECT.nodes(J,2)),N_R(OP.R,d0,abs(NODE.y(ELEMENT_RECT.nodes(J,2))-0.5*PPP),lambda));
                    [X4 Y4 Z4] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODE.x(ELEMENT_RECT.nodes(J,4)),NODE.y(ELEMENT_RECT.nodes(J,4)),N_R(OP.R,d0,abs(NODE.y(ELEMENT_RECT.nodes(J,4))-0.5*PPP),lambda));
                    OP.D1 = sqrt( (X2 - X1)^2 + (Y2 - Y1)^2 + (Z2 - Z1)^2 );
                    OP.D2 = sqrt( (X4 - X1)^2 + (Y4 - Y1)^2 + (Z4 - Z1)^2 );
                end   
                
                MN = mean(NODE.case(ELEMENT_RECT.nodes(J,1:4)));  
                if abs(MN - 2) > 0.0; DEL.E(end+1) = E; end
                if (MN - 2 == 0.0) && (Zn == BaseZ) && (isempty(find(master == rank,1)) == 0); DEL.E(end+1) = E; end
                if (MN > 1) && (MN < 2); ELEMENT_RECT.case(J) = -1; end
                if (MN > 2) && (MN < 3); ELEMENT_RECT.case(J) = 1; end
                if (abs(MN - 2) <= tolS) && (first == 1) && (Zn == BaseZ); first = 0; P.Lefts.ID(end+1) = ELEMENT_RECT.nodes(J,1); end
                if (MN - 2 > 0.0) && (first == 0) && (Zn == BaseZ); P.Rights.ID(end+1) = ELEMENT_RECT.nodes(J,1); first = -1; end
                if (abs(MN - 2) <= tolS) && (first == 1) && (Zn == TopZ1); first = 0; P.Lefts.ID(end+1) = ELEMENT_RECT.nodes(J,4); end
                if (MN - 2 > 0.0) && (first == 0) && (Zn == TopZ1); P.Rights.ID(end+1) = ELEMENT_RECT.nodes(J,4); first = -1; end
                if (ELEMENT_RECT.case(J) == -1) && (Zn <= TopZ1 - tolS); P.Lefts.ID(end+1) = ELEMENT_RECT.nodes(J,3); end
                if (ELEMENT_RECT.case(J) == 1) && (Zn >= BaseZ + tolS); P.Rights.ID(end+1) = ELEMENT_RECT.nodes(J,1); end
            end
        end
    
    elseif GEO.helix.element == 'g' || GEO.helix.element == 'h' || GEO.helix.element == 'j' % If it is an 8-node or 9-node element
        for Zn = BaseZ:2:TopZ2 % For every row but one of vertical nodes
            for Tn = 1:2:RSL.Ttot; E = E + 1; J = J + 1; % For every column but one of circumferential nodes
                if Tn == 1; first = 1; end
                ELEMENT_RECT.index(J) = E; ELEMENT_RECT.type(J) = GEO.helix.element; ELEMENT_RECT.case(J) = 0;
                ELEMENT_RECT.nodes(J,1) = Tn + (Zn-1)*RTcurrent;
                ELEMENT_RECT.nodes(J,2) = Tn + 2 + (Zn-1)*RTcurrent;
                ELEMENT_RECT.nodes(J,3) = Tn + 2 + (Zn+1)*RTcurrent;
                ELEMENT_RECT.nodes(J,4) = Tn + (Zn+1)*RTcurrent;
                ELEMENT_RECT.nodes(J,5) = Tn + 1 + (Zn-1)*RTcurrent;
                ELEMENT_RECT.nodes(J,6) = Tn + 2 + Zn*RTcurrent;
                ELEMENT_RECT.nodes(J,7) = Tn + 1 + (Zn+1)*RTcurrent;
                ELEMENT_RECT.nodes(J,8) = Tn + Zn*RTcurrent;
                if GEO.helix.element == 'j'; ELEMENT_RECT.nodes(J,9) = Tn + 1 + Zn*RTcurrent; end

                if Tn == 1
                    dx = abs(NODE.x(ELEMENT_RECT.nodes(J,1)) - NODE.x(ELEMENT_RECT.nodes(J,2)));
                    dy = abs(NODE.y(ELEMENT_RECT.nodes(J,1)) - NODE.y(ELEMENT_RECT.nodes(J,4)));
                    [X1, Y1, Z1] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODE.x(ELEMENT_RECT.nodes(J,1)),NODE.y(ELEMENT_RECT.nodes(J,1)),N_R(OP.R,d0,abs(NODE.y(ELEMENT_RECT.nodes(J,1))-0.5*PPP),lambda));
                    [X2, Y2, Z2] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODE.x(ELEMENT_RECT.nodes(J,2)),NODE.y(ELEMENT_RECT.nodes(J,2)),N_R(OP.R,d0,abs(NODE.y(ELEMENT_RECT.nodes(J,2))-0.5*PPP),lambda));
                    [X4, Y4, Z4] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODE.x(ELEMENT_RECT.nodes(J,4)),NODE.y(ELEMENT_RECT.nodes(J,4)),N_R(OP.R,d0,abs(NODE.y(ELEMENT_RECT.nodes(J,4))-0.5*PPP),lambda));
                    OP.D1 = sqrt( (X2 - X1)^2 + (Y2 - Y1)^2 + (Z2 - Z1)^2 );
                    OP.D2 = sqrt( (X4 - X1)^2 + (Y4 - Y1)^2 + (Z4 - Z1)^2 );
                end   
                
                MN = mean(NODE.case(ELEMENT_RECT.nodes(J,1:4))); 
                if abs(MN - 2) > 0.0; DEL.E(end+1) = E; end
                if (MN - 2 == 0.0) && (Zn == BaseZ) && (isempty(find(master == rank,1)) == 0); DEL.E(end+1) = E; end
                if (MN > 1) && (MN < 2); ELEMENT_RECT.case(J) = -1; end
                if (MN > 2) && (MN < 3); ELEMENT_RECT.case(J) = 1; end
                if (abs(MN - 2) == 0.0) && (first == 1) && (Zn == BaseZ); first = 0; P.Lefts.ID(end+1) = ELEMENT_RECT.nodes(J,1); end
                if (MN - 2 > 0.0) && (first == 0) && (Zn == BaseZ); P.Rights.ID(end+1) = ELEMENT_RECT.nodes(J,1); first = -1; end
                if (abs(MN - 2) == 0.0) && (first == 1) && (Zn == TopZ2); first = 0; P.Lefts.ID(end+1) = ELEMENT_RECT.nodes(J,4); end
                if (MN - 2 > 0.0) && (first == 0) && (Zn == TopZ2); P.Rights.ID(end+1) = ELEMENT_RECT.nodes(J,4); first = -1; end
                if (ELEMENT_RECT.case(J) == -1) && (Zn < TopZ2); P.Lefts.ID(end+1) = ELEMENT_RECT.nodes(J,3); end
                if (ELEMENT_RECT.case(J) == 1) && (Zn > BaseZ); P.Rights.ID(end+1) = ELEMENT_RECT.nodes(J,1); end
            end
        end
    end

    %%% CALIBRATION OF ELEMENTS
    if isempty(find(slave == rank,1)) == 0; MPI_Send(master,4+G+Gmax+rank,comm,ELEMENT_RECT,P,DEL); end
    if isempty(find(master == rank,1)) == 0 && CPUs > 1
        for S = 1:length(slave)
            [SlaveElements, SlaveP, SlaveDEL] = MPI_Recv(slave(S),4+G+Gmax+slave(S),comm);
            ELEMENT_RECT.index(SlaveElements.index) = SlaveElements.index;
            for L = 1:size(SlaveElements.nodes,2)
                ELEMENT_RECT.nodes(SlaveElements.index,L) = SlaveElements.nodes(:,L);
            end
            ELEMENT_RECT.type(SlaveElements.index) = SlaveElements.type;
            P.Lefts.ID(end+1:end+length(SlaveP.Lefts.ID)) = SlaveP.Lefts.ID;
            P.Lefts2.ID(end+1:end+length(SlaveP.Lefts2.ID)) = SlaveP.Lefts2.ID;
            P.Rights.ID(end+1:end+length(SlaveP.Rights.ID)) = SlaveP.Rights.ID;
            P.Rights2.ID(end+1:end+length(SlaveP.Rights2.ID)) = SlaveP.Rights2.ID;
            P.Ls.x(end+1:end+length(SlaveP.Ls.x)) = SlaveP.Ls.x;
            P.Ls.y(end+1:end+length(SlaveP.Ls.y)) = SlaveP.Ls.y;
            P.Rs.x(end+1:end+length(SlaveP.Rs.x)) = SlaveP.Rs.x;
            P.Rs.y(end+1:end+length(SlaveP.Rs.y)) = SlaveP.Rs.y;
            DEL.N(end+1:end+length(SlaveDEL.N)) = SlaveDEL.N; 
            DEL.E(end+1:end+length(SlaveDEL.E)) = SlaveDEL.E;
        end
    end
    if isempty(find(master == rank,1)) == 0
        LHS = NODE.x(max(P.Lefts.ID)) - P.LTop.x; RHS = P.RTop.x - NODE.x(max(P.Rights.ID));    
        
        if GEO.helix.order == 1
            diff = abs((P.RBot.x - RHS) - (P.LBot.x + LHS))/(length([max(P.Lefts.ID):GEO.helix.order:max(P.Rights.ID)])-1);
            P.Bots.x = [(P.LBot.x + LHS):diff:(P.RBot.x - RHS)]; P.Bots.y = zeros(1,length(P.Bots.x));
            P.Bots.ID = max(NODE.index) + [1:length(P.Bots.x)]; P.Bots.ID(1) = min(P.Lefts.ID); P.Bots.ID(end) = min(P.Rights.ID);      
            DEL.N(end+1:end+length([min(P.Lefts.ID)+1:1:min(P.Rights.ID)-1])) = [min(P.Lefts.ID)+1:1:min(P.Rights.ID)-1];
            for B = 2:length(P.Bots.x)-1
                NODE.index(end+1) = P.Bots.ID(B);
                NODE.x(end+1) = P.Bots.x(B);
                NODE.y(end+1) = P.Bots.y(B);
                NR = N_R(OP.R,d0,abs(P.Bots.y(B)-0.5*PPP),lambda);
                [NODE.X(end+1), NODE.Y(end+1), NODE.Z(end+1)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,P.Bots.x(B),P.Bots.y(B),NR);
            end
            
            P.Tops.ID = [max(P.Lefts.ID):1:max(P.Rights.ID)]; P.Tops.x = NODE.x(P.Tops.ID); P.Tops.y = NODE.y(P.Tops.ID);

            % 1st possible modification
            if length(P.Bots.ID) < length(P.Tops.ID)
                P.Lefts.ID(find(P.Lefts.ID == min(P.Bots.ID))) = min(P.Lefts.ID) - 1;
                P.Lefts.ID(find(P.Lefts.ID == min(P.Bots.ID) + RSL.Ttot + 1)) = [];
                P.Bots.ID = [min(P.Lefts.ID):1:min(P.Rights.ID)]; DEL.N(find(DEL.N == min(P.Bots.ID))) = [];
            end

            P.Ints.ID = [min(P.Lefts.ID):1:min(P.Rights.ID)] + RSL.Ttot + 1; P.Ints.x = NODE.x(P.Ints.ID); P.Ints.y = NODE.y(P.Ints.ID);
            [P.Ints.ID, W] = sort(P.Ints.ID,'descend'); P.Ints.x = P.Ints.x(W); P.Ints.y = P.Ints.y(W);
            
            % Lefts
            Lefts = sort(P.Lefts.ID); P.Lefts.ID = []; P.Lefts.ID(1) = min(Lefts); rem = [];
            for L = 2:length(Lefts) - 1
                if Lefts(L+1) - Lefts(L) == 1; rem(end+1) = Lefts(L);
                else rem(end+1) = Lefts(L); rem = sort(rem,'descend'); P.Lefts.ID(end+1:end+length(rem)) = rem; rem = [];
                end
            end
            P.Lefts.ID(end+1) = max(Lefts); clear('Lefts');
            P.Lefts.x = NODE.x(P.Lefts.ID); P.Lefts.x(1) = P.LBot.x + LHS; P.Lefts.y = NODE.y(P.Lefts.ID);
            Lx = []; Ly = []; LID = []; skip = 0;
            for L = 2:length(P.Lefts.x) 
                if skip == 0; Lx(end+1) = P.Lefts.x(L-1); Ly(end+1) = P.Lefts.y(L-1); LID(end+1) = P.Lefts.ID(L-1); end
                if skip == 1; skip = 0; end
                if abs(P.Lefts.y(L) - P.Lefts.y(L-1)) <= tolS
                    dL = abs(P.Lefts.x(L) - P.Lefts.x(L-1))/dx; n = round(dL);
                    if (abs(dL - n) <= tolS) && (n > 1)
                       if P.Lefts.x(L-1) < P.Lefts.x(L)
                           skip = 1; 
                           Lx(end) = P.Lefts.x(L);
                           Ly(end) = P.Lefts.y(L);
                           LID(end) = P.Lefts.ID(L);
                       end
                       for N = 1:n-1
                           Lx(end+1) = (N/n)*P.Lefts.x(L) + ((n-N)/n)*P.Lefts.x(L-1);
                           Ly(end+1) = P.Lefts.y(L);
                           LID(end+1) = round((N/n)*P.Lefts.ID(L) + ((n-N)/n)*P.Lefts.ID(L-1));
                       end                   
                       if skip == 1
                            Lx(end+1) = P.Lefts.x(L-1);
                            Ly(end+1) = P.Lefts.y(L-1);
                            LID(end+1) = P.Lefts.ID(L-1);
                       end
                    end
                end
                if L == length(P.Lefts.x) && skip == 0; Lx(end+1) = P.Lefts.x(L); Ly(end+1) = P.Lefts.y(L); LID(end+1) = P.Lefts.ID(L); end
            end
            P.Lefts.x = Lx; P.Lefts.y = Ly; P.Lefts.ID = LID;
            
            % Rights
            Rights = sort(P.Rights.ID); P.Rights.ID = []; P.Rights.ID(1) = min(Rights); rem = [];
            for R = 2:length(Rights) - 1
                if Rights(R+1) - Rights(R) == 1; rem(end+1) = Rights(R);
                else rem(end+1) = Rights(R); rem = sort(rem,'descend'); P.Rights.ID(end+1:end+length(rem)) = rem; rem = [];
                end
            end
            P.Rights.ID(end+1) = max(Rights); clear('Rights');
            P.Rights.x = NODE.x(P.Rights.ID); P.Rights.x(1) = P.RBot.x - RHS; P.Rights.y = NODE.y(P.Rights.ID);
            Rx = []; Ry = []; RID = []; skip = 0;
            for R = 2:length(P.Rights.x) 
                if skip == 0; Rx(end+1) = P.Rights.x(R-1); Ry(end+1) = P.Rights.y(R-1); RID(end+1) = P.Rights.ID(R-1); end
                if skip == 1; skip = 0; end
                if abs(P.Rights.y(R) - P.Rights.y(R-1)) < tolS
                    dR = abs(P.Rights.x(R) - P.Rights.x(R-1))/dx; n = round(dR);
                    if (abs(dR - n) <= tolS) && (n > 1)
                       if P.Rights.x(R-1) < P.Rights.x(R)
                           skip = 1; 
                           Rx(end) = P.Rights.x(R);
                           Ry(end) = P.Rights.y(R);
                           RID(end) = P.Rights.ID(R);
                       end
                       for N = 1:n-1
                           Rx(end+1) = (N/n)*P.Rights.x(R) + ((n-N)/n)*P.Rights.x(R-1);
                           Ry(end+1) = P.Rights.y(R);
                           RID(end+1) = round((N/n)*P.Rights.ID(R) + ((n-N)/n)*P.Rights.ID(R-1));
                       end                    
                       if skip == 1
                            Rx(end+1) = P.Rights.x(R-1);
                            Ry(end+1) = P.Rights.y(R-1);
                            RID(end+1) = P.Rights.ID(R-1);
                       end
                    end
                end
                if R == length(P.Rights.x) && skip == 0; Rx(end+1) = P.Rights.x(R); Ry(end+1) = P.Rights.y(R); RID(end+1) = P.Rights.ID(R); end
            end
            P.Rights.x = Rx; P.Rights.y = Ry; P.Rights.ID = RID;
            
        elseif GEO.helix.order == 2
            % Centres
            P.Ints.ID = [min(P.Lefts.ID):2:min(P.Rights.ID)] + 2*RSL.Ttot + 2;
            P.Ints.x = NODE.x(P.Ints.ID); P.Ints.y = NODE.y(P.Ints.ID);
            [P.Ints.ID, W] = sort(P.Ints.ID,'descend'); P.Ints.x = P.Ints.x(W); P.Ints.y = P.Ints.y(W);
            
            P.Ints2.ID = [min(P.Lefts.ID)+1:2:min(P.Rights.ID)-1] + 2*RSL.Ttot + 2;
            P.Ints2.x = NODE.x(P.Ints2.ID); P.Ints2.y = NODE.y(P.Ints2.ID);
            [P.Ints2.ID, W] = sort(P.Ints2.ID,'descend'); P.Ints2.x = P.Ints2.x(W); P.Ints2.y = P.Ints2.y(W);
            
            P.Tops.ID = [max(P.Lefts.ID):2:max(P.Rights.ID)]; 
            P.Tops.x = NODE.x(P.Tops.ID); P.Tops.y = NODE.y(P.Tops.ID);
            
            P.Tops2.ID = [max(P.Lefts.ID)+1:2:max(P.Rights.ID)-1];
            P.Tops2.x = NODE.x(P.Tops2.ID); P.Tops2.y = NODE.y(P.Tops2.ID);     
            
            P.Bots.x = ((P.Tops.x-min(P.Tops.x))+(P.LBot.x + LHS)); P.Bots.y = zeros(1,length(P.Bots.x));            
            P.Bots.ID = max(NODE.index) + [1:length(P.Bots.x)]; P.Bots.ID(1) = min(P.Lefts.ID); P.Bots.ID(end) = min(P.Rights.ID);
            DEL.N(end+1:end+length([min(P.Lefts.ID)+2:2:min(P.Rights.ID)-2])) = [min(P.Lefts.ID)+2:2:min(P.Rights.ID)-2];
            for B = 2:length(P.Bots.x)-1
                NODE.index(end+1) = P.Bots.ID(B);
                NODE.x(end+1) = P.Bots.x(B);
                NODE.y(end+1) = P.Bots.y(B);
                NR = N_R(OP.R,d0,abs(P.Bots.y(B)-0.5*PPP),lambda);
                [NODE.X(end+1), NODE.Y(end+1), NODE.Z(end+1)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,P.Bots.x(B),P.Bots.y(B),NR);
            end
            DEL.N(end+1:end+length([min(P.Lefts.ID)+2:2:min(P.Rights.ID)-2])) = [min(P.Lefts.ID)+2:2:min(P.Rights.ID)-2]+ RSL.Ttot +1;
            
            P.Bots2.x = ((P.Tops2.x-min(P.Tops2.x))+(P.LBot.x + LHS + abs(min(P.Tops2.x) - min(P.Tops.x)))); P.Bots2.y = zeros(1,length(P.Bots2.x));
            P.Bots2.ID = max(P.Bots.ID) + [1:length(P.Bots2.x)];
            DEL.N(end+1:end+length([min(P.Lefts.ID)+1:2:min(P.Rights.ID)-1])) = [min(P.Lefts.ID)+1:2:min(P.Rights.ID)-1];
            for B2 = 1:length(P.Bots2.x)
                NODE.index(end+1) = P.Bots2.ID(B2);
                NODE.x(end+1) = P.Bots2.x(B2);
                NODE.y(end+1) = P.Bots2.y(B2);
                NR = N_R(OP.R,d0,abs(P.Bots2.y(B2)-0.5*PPP),lambda);
                [NODE.X(end+1), NODE.Y(end+1), NODE.Z(end+1)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,P.Bots2.x(B2),P.Bots2.y(B2),NR);
            end
            DEL.N(end+1:end+length([min(P.Lefts.ID)+1:2:min(P.Rights.ID)-1])) = [min(P.Lefts.ID)+1:2:min(P.Rights.ID)-1]+ RSL.Ttot +1;
            
            % Lefts
            Lefts = sort(P.Lefts.ID); P.Lefts.ID = []; P.Lefts.ID(1) = min(Lefts); rem = [];
            for L = 2:length(Lefts) - 1
                if Lefts(L+1) - Lefts(L) == 2; rem(end+1) = Lefts(L);
                else rem(end+1) = Lefts(L); rem = sort(rem,'descend'); P.Lefts.ID(end+1:end+length(rem)) = rem; rem = [];
                end
            end
            P.Lefts.ID(end+1) = max(Lefts); clear('Lefts');
            P.Lefts.x = NODE.x(P.Lefts.ID); P.Lefts.y = NODE.y(P.Lefts.ID);
            NR = N_R(OP.R,d0,abs(NODE.y(P.Lefts.ID(1))-0.5*PPP),lambda);
            P.Lefts.x(1) = P.LBot.x + LHS; [NODE.X(P.Lefts.ID(1)), NODE.Y(P.Lefts.ID(1)), NODE.Z(P.Lefts.ID(1))] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,P.Lefts.x(1),P.Lefts.y(1),NR);                      
            Lx = []; Ly = []; LID = []; skip = 0;
            for L = 2:length(P.Lefts.x)
                if skip == 0; Lx(end+1) = P.Lefts.x(L-1); Ly(end+1) = P.Lefts.y(L-1); LID(end+1) = P.Lefts.ID(L-1); end
                if skip == 1; skip = 0; end
                if abs(P.Lefts.y(L) - P.Lefts.y(L-1)) <= tolS
                    dL = abs(P.Lefts.x(L) - P.Lefts.x(L-1))/dx; n = round(dL);
                    if (abs(dL - n) <= tolS) && (n > 1)
                       if P.Lefts.x(L-1) < P.Lefts.x(L)
                           skip = 1; 
                           Lx(end) = P.Lefts.x(L);
                           Ly(end) = P.Lefts.y(L);
                           LID(end) = P.Lefts.ID(L);
                       end
                       for N = 1:n-1
                           Lx(end+1) = (N/n)*P.Lefts.x(L) + ((n-N)/n)*P.Lefts.x(L-1);
                           Ly(end+1) = P.Lefts.y(L);
                           LID(end+1) = round((N/n)*P.Lefts.ID(L) + ((n-N)/n)*P.Lefts.ID(L-1));
                       end                   
                       if skip == 1
                            Lx(end+1) = P.Lefts.x(L-1);
                            Ly(end+1) = P.Lefts.y(L-1);
                            LID(end+1) = P.Lefts.ID(L-1);
                       end
                    end
                end
                if L == length(P.Lefts.x) && skip == 0; Lx(end+1) = P.Lefts.x(L); Ly(end+1) = P.Lefts.y(L); LID(end+1) = P.Lefts.ID(L); end
            end
            P.Lefts.x = Lx; P.Lefts.y = Ly; P.Lefts.ID = LID;           

            for L = 1:length(P.Lefts.ID) - 1
                P.Lefts2.ID(end+1) = round((P.Lefts.ID(L) + P.Lefts.ID(L+1))/2);
                P.Lefts2.x(end+1) = (P.Lefts.x(L) + P.Lefts.x(L+1))/2;
                P.Lefts2.y(end+1) = (P.Lefts.y(L) + P.Lefts.y(L+1))/2;
                NR = N_R(OP.R,d0,abs(NODE.y(P.Lefts2.ID(end))-0.5*PPP),lambda);
                [NODE.X(P.Lefts2.ID(end)), NODE.Y(P.Lefts2.ID(end)), NODE.Z(P.Lefts2.ID(end))] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,P.Lefts2.x(end),P.Lefts2.y(end),NR);
            end
            NODE.x(P.Lefts.ID) = P.Lefts.x; NODE.x(P.Lefts2.ID) = P.Lefts2.x;
            for L = 1:length(P.Lefts.ID)
                NR = N_R(OP.R,d0,abs(NODE.y(P.Lefts.ID(L))-0.5*PPP),lambda);
                [NODE.X(P.Lefts.ID(L)), NODE.Y(P.Lefts.ID(L)), NODE.Z(P.Lefts.ID(L))] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,P.Lefts.x(L),P.Lefts.y(L),NR);
            end

            % Rights
            Rights = sort(P.Rights.ID); P.Rights.ID = []; P.Rights.ID(1) = min(Rights); rem = [];
            for R = 2:length(Rights)-1
                if Rights(R+1) - Rights(R) == 2; rem(end+1) = Rights(R);
                else rem(end+1) = Rights(R); rem = sort(rem,'descend'); P.Rights.ID(end+1:end+length(rem)) = rem; rem = [];
                end
            end
            P.Rights.ID(end+1) = max(Rights); clear('Rights');
            P.Rights.x = NODE.x(P.Rights.ID); P.Rights.y = NODE.y(P.Rights.ID);
            NR = N_R(OP.R,d0,abs(NODE.y(P.Rights.ID(1))-0.5*PPP),lambda);
            P.Rights.x(1) = P.RBot.x - RHS; [NODE.X(P.Rights.ID(1)), NODE.Y(P.Rights.ID(1)), NODE.Z(P.Rights.ID(1))] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,P.Rights.x(1),P.Rights.y(1),NR);
           Rx = []; Ry = []; RID = []; skip = 0;
            for R = 2:length(P.Rights.x)
                if skip == 0; Rx(end+1) = P.Rights.x(R-1); Ry(end+1) = P.Rights.y(R-1); RID(end+1) = P.Rights.ID(R-1); end
                if skip == 1; skip = 0; end
                if abs(P.Rights.y(R) - P.Rights.y(R-1)) < tolS
                    dR = abs(P.Rights.x(R) - P.Rights.x(R-1))/dx; n = round(dR);
                    if (abs(dR - n) <= tolS) && (n > 1)
                       if P.Rights.x(R-1) < P.Rights.x(R)
                           skip = 1; 
                           Rx(end) = P.Rights.x(R);
                           Ry(end) = P.Rights.y(R);
                           RID(end) = P.Rights.ID(R);
                       end
                       for N = 1:n-1
                           Rx(end+1) = (N/n)*P.Rights.x(R) + ((n-N)/n)*P.Rights.x(R-1);
                           Ry(end+1) = P.Rights.y(R);
                           RID(end+1) = round((N/n)*P.Rights.ID(R) + ((n-N)/n)*P.Rights.ID(R-1));
                       end                    
                       if skip == 1
                            Rx(end+1) = P.Rights.x(R-1);
                            Ry(end+1) = P.Rights.y(R-1);
                            RID(end+1) = P.Rights.ID(R-1);
                       end
                    end
                end
                if R == length(P.Rights.x) && skip == 0; Rx(end+1) = P.Rights.x(R); Ry(end+1) = P.Rights.y(R); RID(end+1) = P.Rights.ID(R); end
            end
            P.Rights.x = Rx; P.Rights.y = Ry; P.Rights.ID = RID;
            
            for R = 1:length(P.Rights.ID) - 1
                P.Rights2.ID(end+1) = (P.Rights.ID(R) + P.Rights.ID(R+1))/2;
                P.Rights2.x(end+1) = (P.Rights.x(R) + P.Rights.x(R+1))/2;
                P.Rights2.y(end+1) = (P.Rights.y(R) + P.Rights.y(R+1))/2;
                NR = N_R(OP.R,d0,abs(NODE.y(P.Rights2.ID(end))-0.5*PPP),lambda);
                [NODE.X(P.Rights2.ID(end)), NODE.Y(P.Rights2.ID(end)), NODE.Z(P.Rights2.ID(end))] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,P.Rights2.x(end),P.Rights2.y(end),NR);
            end
            NODE.x(P.Rights.ID) = P.Rights.x; NODE.x(P.Rights2.ID) = P.Rights2.x;
            for R = 1:length(P.Rights.ID)
                NR = N_R(OP.R,d0,abs(NODE.y(P.Rights.ID(R))-0.5*PPP),lambda);
                [NODE.X(P.Rights.ID(R)), NODE.Y(P.Rights.ID(R)), NODE.Z(P.Rights.ID(R))] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,P.Rights.x(R),P.Rights.y(R),NR);
            end
        end
        P.LTop.ID = max(NODE.index)+1; P.LBot.ID = P.LTop.ID + 1;
        P.RTop.ID = P.LTop.ID*(1+CPUs); P.RBot.ID = P.RTop.ID + 1;
    end
    OP.P = P; OP.DEL = DEL;
    ELEMENT_RECT = rmfield(ELEMENT_RECT,'case'); NODE = rmfield(NODE,'case');
end