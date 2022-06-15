function [NODES,ELEMENTS,OP,NSET,ELSET] = AJSABQ_triangulation(PP,RSL,OP,NODES,ELEMENTS,GEO,MATERIAL)
% Triangulation preparation algorithm
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:27 (previously 06/09/12 - 15:08)

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

master = PP.master; slave = PP.slave; rank = PP.rank; comm = PP.comm; CPUs = PP.CPUs; conformal = GEO.helix.conformal;
Pitch = GEO.helix.pitch; np = GEO.helix.no_p; alpha = GEO.helix.incline*pi/180; Alpha = GEO.helix.incline;

% Preliminaries
fL_x =  @(y,PPP,m)  (PPP - y)/m;
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

if isempty(find(master == rank,1)) == 0
    no = OP.Tols.sideE;
    xoff = no*abs(RSL.t2(1) - RSL.t1(1))/RSL.T(1); OP.DEL.E = unique(OP.DEL.E); OP.DEL.N = unique(OP.DEL.N);
    ELEMENTS.index(OP.DEL.E) = []; ELEMENTS.type(OP.DEL.E) = []; ELEMENTS.nodes(OP.DEL.E,:) = []; OP.Hel_Quad = length(ELEMENTS.index);
    NODES.index(OP.DEL.N) = []; NODES.x(OP.DEL.N) = []; NODES.y(OP.DEL.N) = [];
    NODES.X(OP.DEL.N) = []; NODES.Y(OP.DEL.N) = []; NODES.Z(OP.DEL.N) = [];

    mn = find(NODES.x == min(NODES.x)); if length(mn) == 1; NODES.x(mn) = []; NODES.y(mn) = []; NODES.X(mn) = []; NODES.Y(mn) = []; NODES.Z(mn) = []; NODES.index(mn) = []; end
    mx = find(NODES.x == max(NODES.x)); if length(mx) == 1; NODES.x(mx) = []; NODES.y(mx) = []; NODES.X(mx) = []; NODES.Y(mx) = []; NODES.Z(mx) = []; NODES.index(mx) = []; end

    NSET = struct('Left',[],'Top',[],'Right',[],'Bottom',[]);
    ELSET = struct('Left',[],'Top',[],'Right',[],'Bottom',[]);
end

if CPUs > 3 && rank > 2; ELSET = []; NSET = []; end

switch CPUs
    case 1
        if isempty(find(master == rank,1)) == 0
            % CPU = 1 - Master Thread - Triangulating the LHS strip
            BL_x = []; BL_y = []; BL_ID = [];          
            BL_x(1) = OP.P.LTop.x; BL_y(1) = OP.P.LTop.y; BL_ID(1) = max(NODES.index)+1;
            N = round(2*pi*OP.R/OP.D1);
            for Y = 1:N-1
                BL_x(end+1) = OP.P.LTop.x + Y*(OP.P.LBot.x - OP.P.LTop.x)/N;
                BL_y(end+1) = OP.P.LTop.y + Y*(OP.P.LBot.y - OP.P.LTop.y)/N;
                if OP.Tols.interior == 1
                    OP.P.Ls.x(end+1) = fL_x(BL_y(end),PPP,m) + (1/(OP.Tols.sideE+1))*OP.xoff;
                    OP.P.Ls.y(end+1) = BL_y(end);
                end
            end
            BL_ID(1:length(BL_x)) = max(BL_ID) + [1:length(BL_x)];
            BL_x(end+1) = OP.P.LBot.x; BL_y(end+1) = OP.P.LBot.y; BL_ID(end+1) = max(BL_ID) + 1;
            BL_x(end+1) = 0.5*(min(OP.P.Bots.x) + OP.P.LBot.x); BL_y(end+1) = OP.P.LBot.y; BL_ID(end+1) = max(BL_ID) + 1; bts(1:3) = [BL_ID(end-1:end) OP.P.Bots.ID(1)];
            BL_x(end+1:end+length(OP.P.Lefts.ID)) = OP.P.Lefts.x; BL_y(end+1:end+length(OP.P.Lefts.ID)) = OP.P.Lefts.y; BL_ID(end+1:end+length(OP.P.Lefts.ID)) = OP.P.Lefts.ID;
            BL_x(end+1) = 0.5*(min(OP.P.Tops.x) + OP.P.LTop.x); BL_y(end+1) = OP.P.LTop.y; BL_ID(end+1) = max(BL_ID) + 1; tps(1:3) = [BL_ID(1) BL_ID(end) OP.P.Tops.ID(1)];
            if OP.Tols.interior == 1
                IL_x = OP.P.Ls.x; IL_y = OP.P.Ls.y; IL_ID = [];
                IL_ID(1:length(IL_x)) = max(BL_ID) + [1:length(IL_x)];
            else IL_x = []; IL_y = []; IL_ID = []; end

            [TRI_L,VER_L,TRI_L_boundary,VER_L_boundary] = AJSABQ_Delaunay(BL_x,BL_y,BL_ID,IL_x,IL_y,IL_ID,OP.Tols.tolRND,OP.Tols.tolNN,OP.Tols.tolBG,GEO.helix.order,max(NODES.index)+1,max(ELEMENTS.index)+1,OP.Tols.scale);
            disp('Consolidating');
            if GEO.helix.order == 2
                ns = []; bt = find(VER_L.y < min(VER_L.y) + OP.Tols.tolNN); tp = find(VER_L.y > max(VER_L.y) - OP.Tols.tolNN); 
                bts(end+1:end+length(bt)) = VER_L.ID(bt); tps(end+1:end+length(tp)) = VER_L.ID(tp);
                bts = unique(bts); tps = unique(tps); bts_x = []; tps_x = []; 
                for bt = 1:length(bts)
                    bts_x(bt) = VER_L.x(find(VER_L.ID == bts(bt)));
                    tps_x(bt) = VER_L.x(find(VER_L.ID == tps(bt)));
                end
                [bts_x, b] = sort(bts_x); [tps_x, t] = sort(tps_x);
                bts = bts(b); tps = tps(t);
            end
            for T = 1:length(TRI_L.ID)
                n1 = find(VER_L.ID == TRI_L.vertices(T,1)); x1 = VER_L.x(n1); y1 = VER_L.y(n1);
                n2 = find(VER_L.ID == TRI_L.vertices(T,2)); x2 = VER_L.x(n2); y2 = VER_L.y(n2);
                n3 = find(VER_L.ID == TRI_L.vertices(T,3)); x3 = VER_L.x(n3); y3 = VER_L.y(n3);
                if GEO.helix.order == 2; WH = [];
                    n4 = find(VER_L.ID == TRI_L.vertices(T,4)); x4 = VER_L.x(n4); y4 = VER_L.y(n4); b4 = 1; if isempty(find(VER_L_boundary == VER_L.ID(n4), 1)); b4 = 0; end
                    n5 = find(VER_L.ID == TRI_L.vertices(T,5)); x5 = VER_L.x(n5); y5 = VER_L.y(n5); b5 = 1; if isempty(find(VER_L_boundary == VER_L.ID(n5), 1)); b5 = 0; end
                    n6 = find(VER_L.ID == TRI_L.vertices(T,6)); x6 = VER_L.x(n6); y6 = VER_L.y(n6); b6 = 1; if isempty(find(VER_L_boundary == VER_L.ID(n6), 1)); b6 = 0; end
                    if b4 == 1; WH(end+1) = 12; end
                    if b5 == 1; WH(end+1) = 23; end
                    if b6 == 1; WH(end+1) = 13; end
                    for W = 1:length(WH)
                        if WH(W) == 12; x_e = x4; y_e = y4; n = n4; N = 4; end
                        if WH(W) == 13; x_e = x6; y_e = y6; n = n6; N = 6; end
                        if WH(W) == 23; x_e = x5; y_e = y5; n = n5; N = 5; end
                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Lefts2.x; OP.P.Lefts2.y]);
                        PTx = OP.P.Lefts2.x(IDX); PTy = OP.P.Lefts2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_L.vertices(T,N) = OP.P.Lefts2.ID(IDX); ns(end+1) = n; end
                    end
                end
                x_edge12 = 0.5*(x1 + x2); y_edge12 = 0.5*(y1 + y2);
                x_edge13 = 0.5*(x1 + x3); y_edge13 = 0.5*(y1 + y3);
                x_edge23 = 0.5*(x2 + x3); y_edge23 = 0.5*(y2 + y3);
                if abs(x_edge12 - fL_x(y_edge12,PPP,m)) < OP.Tols.tolNN
                    NSET.Left(end+1) = VER_L.ID(n1); NSET.Left(end+1) = VER_L.ID(n2); ELSET.Left(end+1) = TRI_L.ID(T);
                    if GEO.helix.order == 1; TRI_L.vertices(T,:) = [VER_L.ID(n2) VER_L.ID(n3) VER_L.ID(n1)]; end
                    if GEO.helix.order == 2; TRI_L.vertices(T,:) = [VER_L.ID(n2) VER_L.ID(n3) VER_L.ID(n1) VER_L.ID(n5) VER_L.ID(n6) VER_L.ID(n4)]; NSET.Left(end+1) = VER_L.ID(n4); end
                end
                if abs(x_edge13 - fL_x(y_edge13,PPP,m)) < OP.Tols.tolNN
                    NSET.Left(end+1) = VER_L.ID(n1); NSET.Left(end+1) = VER_L.ID(n3); ELSET.Left(end+1) = TRI_L.ID(T);
                    if GEO.helix.order == 2; NSET.Left(end+1) = VER_L.ID(n6); end
                end
                if abs(x_edge23 - fL_x(y_edge23,PPP,m)) < OP.Tols.tolNN
                    NSET.Left(end+1) = VER_L.ID(n2); NSET.Left(end+1) = VER_L.ID(n3); ELSET.Left(end+1) = TRI_L.ID(T);
                    if GEO.helix.order == 1; TRI_L.vertices(T,:) = [VER_L.ID(n3) VER_L.ID(n1) VER_L.ID(n2)]; end
                    if GEO.helix.order == 2; TRI_L.vertices(T,:) = [VER_L.ID(n3) VER_L.ID(n1) VER_L.ID(n2) VER_L.ID(n6) VER_L.ID(n4) VER_L.ID(n5)]; NSET.Left(end+1) = VER_L.ID(n5); end
                end
                for L = 1:length(bts)
                    n = find(TRI_L.vertices(T,:) == bts(L));
                    if ~ isempty(n); TRI_L.vertices(T,n) = tps(L); end
                end
                if (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) < 0.0 % Verifying the element normal
                   if GEO.helix.order == 1; TRI_L.vertices(T,:) = [TRI_L.vertices(T,3) TRI_L.vertices(T,2) TRI_L.vertices(T,1)]; end
                   if GEO.helix.order == 2; TRI_L.vertices(T,:) = [TRI_L.vertices(T,3) TRI_L.vertices(T,2) TRI_L.vertices(T,1) TRI_L.vertices(T,5) TRI_L.vertices(T,4) TRI_L.vertices(T,6)]; end                   
                end
            end
            if GEO.helix.order == 1; en = 1; end
            if GEO.helix.order == 2; VER_L.ID(ns) = []; VER_L.x(ns) = []; VER_L.y(ns) = [];
                if GEO.helix.element == 'j'; en = 3; else en = 2; end
            end
            NSET.Bottom = unique(NSET.Left); NSET = rmfield(NSET,'Left');
            ELSET.Bottom = ELSET.Left; ELSET = rmfield(ELSET,'Left');

            ELEMENTS.nodes(end+1:end+size(TRI_L.vertices,1),:) = [TRI_L.vertices zeros(size(TRI_L.vertices,1),en)];
            ELEMENTS.index(end+1:end+length(TRI_L.ID)) = TRI_L.ID;
            if GEO.helix.order == 1; ELEMENTS.type(end+1:end+length(TRI_L.ID)) = char(ones(1,length(TRI_L.ID))*'d'); len = length(NODES.index); end
            if GEO.helix.order == 2; ELEMENTS.type(end+1:end+length(TRI_L.ID)) = char(ones(1,length(TRI_L.ID))*'f'); len = length(NODES.index); end
            NODES.index(end+1:end+length(VER_L.ID)) = VER_L.ID;
            NODES.x(end+1:end+length(VER_L.ID)) = VER_L.x;
            NODES.y(end+1:end+length(VER_L.ID)) = VER_L.y;
            for N = 1:length(VER_L.ID)
                % Radial coordinate
                NR = N_R(OP.R,d0,abs(NODES.y(len+N)-0.5*PPP),lambda);
                [NODES.X(len+N), NODES.Y(len+N), NODES.Z(len+N)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODES.x(len+N),NODES.y(len+N),NR);
            end


            % CPU = 1 - Master Thread - Triangulating the RHS strip
            BR_x = []; BR_y = []; BR_ID = []; clear('bts','tps');        
            BR_x(1) = OP.P.RTop.x; BR_y(1) = OP.P.RTop.y; BR_ID(1) = max(NODES.index)+1;
            N = round(2*pi*OP.R/OP.D1);
            for Y = 1:N-1
                BR_x(end+1) = OP.P.RTop.x + Y*(OP.P.RBot.x - OP.P.RTop.x)/N;
                BR_y(end+1) = OP.P.RTop.y + Y*(OP.P.RBot.y - OP.P.RTop.y)/N;
                if OP.Tols.interior == 1
                    OP.P.Rs.x(end+1) = fR_x(BR_y(end),x_max,PPP,m) - (1/(OP.Tols.sideE+1))*OP.xoff;
                    OP.P.Rs.y(end+1) = BR_y(end);
                end
            end
            BR_ID(1:length(BR_x)) = max(BR_ID) + [1:length(BR_x)];
            BR_x(end+1) = OP.P.RBot.x; BR_y(end+1) = OP.P.RBot.y; BR_ID(end+1) = max(BR_ID) + 1;
            BR_x(end+1) = 0.5*(max(OP.P.Bots.x) + OP.P.RBot.x); BR_y(end+1) = OP.P.RBot.y; BR_ID(end+1) = max(BR_ID) + 1; bts(1:3) = [BR_ID(end-1:end) OP.P.Bots.ID(end)]; bts = fliplr(bts);
            BR_x(end+1:end+length(OP.P.Rights.ID)) = OP.P.Rights.x; BR_y(end+1:end+length(OP.P.Rights.ID)) = OP.P.Rights.y; BR_ID(end+1:end+length(OP.P.Rights.ID)) = OP.P.Rights.ID;
            BR_x(end+1) = 0.5*(max(OP.P.Tops.x) + OP.P.RTop.x); BR_y(end+1) = OP.P.RTop.y; BR_ID(end+1) = max(BR_ID) + 1; tps(1:3) = [BR_ID(1) BR_ID(end) OP.P.Tops.ID(end)]; tps = fliplr(tps);
            if OP.Tols.interior == 1
                IR_x = OP.P.Rs.x; IR_y = OP.P.Rs.y; IR_ID = [];
                IR_ID(1:length(IR_x)) = max(BR_ID) + [1:length(IR_x)];
            else IR_x = []; IR_y = []; IR_ID = []; end

            [TRI_R,VER_R,TRI_R_boundary,VER_R_boundary] = AJSABQ_Delaunay(BR_x,BR_y,BR_ID,IR_x,IR_y,IR_ID,OP.Tols.tolRND,OP.Tols.tolNN,OP.Tols.tolBG,GEO.helix.order,max(NODES.index)+1,max(ELEMENTS.index)+1,OP.Tols.scale);
            disp('Consolidating');
            if GEO.helix.order == 2
                ns = []; bt = find(VER_R.y < min(VER_R.y) + OP.Tols.tolNN); tp = find(VER_R.y > max(VER_R.y) - OP.Tols.tolNN); 
                bts(end+1:end+length(bt)) = VER_R.ID(bt); tps(end+1:end+length(tp)) = VER_R.ID(tp);
                bts = unique(bts); tps = unique(tps); bts_x = []; tps_x = [];
                for bt = 1:length(bts)
                    bts_x(bt) = VER_R.x(find(VER_R.ID == bts(bt)));
                    tps_x(bt) = VER_R.x(find(VER_R.ID == tps(bt)));
                end
                [bts_x, b] = sort(bts_x); [tps_x, t] = sort(tps_x);
                bts = bts(b); tps = tps(t);
            end
            for T = 1:length(TRI_R.ID)
                n1 = find(VER_R.ID == TRI_R.vertices(T,1)); x1 = VER_R.x(n1); y1 = VER_R.y(n1);
                n2 = find(VER_R.ID == TRI_R.vertices(T,2)); x2 = VER_R.x(n2); y2 = VER_R.y(n2);
                n3 = find(VER_R.ID == TRI_R.vertices(T,3)); x3 = VER_R.x(n3); y3 = VER_R.y(n3);
                if GEO.helix.order == 2; WH = [];
                    n4 = find(VER_R.ID == TRI_R.vertices(T,4)); x4 = VER_R.x(n4); y4 = VER_R.y(n4); b4 = 1; if isempty(find(VER_R_boundary == VER_R.ID(n4), 1)); b4 = 0; end
                    n5 = find(VER_R.ID == TRI_R.vertices(T,5)); x5 = VER_R.x(n5); y5 = VER_R.y(n5); b5 = 1; if isempty(find(VER_R_boundary == VER_R.ID(n5), 1)); b5 = 0; end
                    n6 = find(VER_R.ID == TRI_R.vertices(T,6)); x6 = VER_R.x(n6); y6 = VER_R.y(n6); b6 = 1; if isempty(find(VER_R_boundary == VER_R.ID(n6), 1)); b6 = 0; end
                    if b4 == 1; WH(end+1) = 12; end
                    if b5 == 1; WH(end+1) = 23; end
                    if b6 == 1; WH(end+1) = 13; end
                    for W = 1:length(WH)
                        if WH(W) == 12; x_e = x4; y_e = y4; n = n4; N = 4; end
                        if WH(W) == 13; x_e = x6; y_e = y6; n = n6; N = 6; end
                        if WH(W) == 23; x_e = x5; y_e = y5; n = n5; N = 5; end
                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Rights2.x; OP.P.Rights2.y]);
                        PTx = OP.P.Rights2.x(IDX); PTy = OP.P.Rights2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_R.vertices(T,N) = OP.P.Rights2.ID(IDX); ns(end+1) = n; end
                    end
                end
                x_edge12 = 0.5*(x1 + x2); y_edge12 = 0.5*(y1 + y2);
                x_edge13 = 0.5*(x1 + x3); y_edge13 = 0.5*(y1 + y3);
                x_edge23 = 0.5*(x2 + x3); y_edge23 = 0.5*(y2 + y3);
                if abs(x_edge12 - fR_x(y_edge12,x_max,PPP,m)) < OP.Tols.tolNN
                    NSET.Right(end+1) = VER_R.ID(n1); NSET.Right(end+1) = VER_R.ID(n2); ELSET.Right(end+1) = TRI_R.ID(T);
                    if GEO.helix.order == 1; TRI_R.vertices(T,:) = [VER_R.ID(n2) VER_R.ID(n3) VER_R.ID(n1)]; end
                    if GEO.helix.order == 2; TRI_R.vertices(T,:) = [VER_R.ID(n2) VER_R.ID(n3) VER_R.ID(n1) VER_R.ID(n5) VER_R.ID(n6) VER_R.ID(n4)]; NSET.Right(end+1) = VER_R.ID(n4); end
                end
                if abs(x_edge13 - fR_x(y_edge13,x_max,PPP,m)) < OP.Tols.tolNN
                    NSET.Right(end+1) = VER_R.ID(n1); NSET.Right(end+1) = VER_R.ID(n3); ELSET.Right(end+1) = TRI_R.ID(T);
                    if GEO.helix.order == 2; NSET.Right(end+1) = VER_R.ID(n6); end
                end
                if abs(x_edge23 - fR_x(y_edge23,x_max,PPP,m)) < OP.Tols.tolNN
                    NSET.Right(end+1) = VER_R.ID(n2); NSET.Right(end+1) = VER_R.ID(n3); ELSET.Right(end+1) = TRI_R.ID(T);
                    if GEO.helix.order == 1; TRI_R.vertices(T,:) = [VER_R.ID(n3) VER_R.ID(n1) VER_R.ID(n2)]; end
                    if GEO.helix.order == 2; TRI_R.vertices(T,:) = [VER_R.ID(n3) VER_R.ID(n1) VER_R.ID(n2) VER_R.ID(n6) VER_R.ID(n4) VER_R.ID(n5)]; NSET.Right(end+1) = VER_R.ID(n5); end
                end
                for R = 1:length(bts)
                    n = find(TRI_R.vertices(T,:) == bts(R));
                    if ~ isempty(n); TRI_R.vertices(T,n) = tps(R); end
                end
                if (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) < 0.0 % Verifying the consistency of the element normal
                   if GEO.helix.order == 1; TRI_R.vertices(T,:) = [TRI_R.vertices(T,3) TRI_R.vertices(T,2) TRI_R.vertices(T,1)]; end
                   if GEO.helix.order == 2; TRI_R.vertices(T,:) = [TRI_R.vertices(T,3) TRI_R.vertices(T,2) TRI_R.vertices(T,1) TRI_R.vertices(T,5) TRI_R.vertices(T,4) TRI_R.vertices(T,6)]; end                   
                end
            end
            if GEO.helix.order == 1; en = 1; end
            if GEO.helix.order == 2
                VER_R.ID(ns) = []; VER_R.x(ns) = []; VER_R.y(ns) = [];
                if GEO.helix.element == 'j'; en = 3; else en = 2; end
            end
            NSET.Top = unique(NSET.Right); NSET = rmfield(NSET,'Right');
            ELSET.Top = ELSET.Right; ELSET = rmfield(ELSET,'Right');

            ELEMENTS.nodes(end+1:end+size(TRI_R.vertices,1),:) = [TRI_R.vertices zeros(size(TRI_R.vertices,1),en)];
            ELEMENTS.index(end+1:end+length(TRI_R.ID)) = TRI_R.ID;
            if GEO.helix.order == 1; ELEMENTS.type(end+1:end+length(TRI_R.ID)) = char(ones(1,length(TRI_R.ID))*'d'); len = length(NODES.index); end
            if GEO.helix.order == 2; ELEMENTS.type(end+1:end+length(TRI_R.ID)) = char(ones(1,length(TRI_R.ID))*'f'); len = length(NODES.index); end
            NODES.index(end+1:end+length(VER_R.ID)) = VER_R.ID;
            NODES.x(end+1:end+length(VER_R.ID)) = VER_R.x;
            NODES.y(end+1:end+length(VER_R.ID)) = VER_R.y;
            for N = 1:length(VER_R.ID)
                % Radial coordinate
                NR = N_R(OP.R,d0,abs(NODES.y(len+N)-0.5*PPP),lambda);
                [NODES.X(len+N), NODES.Y(len+N), NODES.Z(len+N)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODES.x(len+N),NODES.y(len+N),NR);
            end


            % CPU = 1 - Master Thread - Triangulating the misalingment
            BB_ID = [OP.P.Bots.ID OP.P.Ints.ID]; BB_x = [OP.P.Bots.x OP.P.Ints.x]; BB_y = [OP.P.Bots.y OP.P.Ints.y]; ns = [];

            [TRI_B,VER_B,TRI_B_boundary,VER_B_boundary] = AJSABQ_Delaunay(BB_x,BB_y,BB_ID,[],[],[],OP.Tols.tolRND,OP.Tols.tolNN,OP.Tols.tolBG,GEO.helix.order,max(NODES.index)+1,max(ELEMENTS.index)+1,0);
            disp('Consolidating');
            for T = 1:length(TRI_B.ID)
                if GEO.helix.order == 2; WH = [];
                    n4 = find(VER_B.ID == TRI_B.vertices(T,4)); x4 = VER_B.x(n4); y4 = VER_B.y(n4); b4 = 1; if isempty(find(VER_B_boundary == VER_B.ID(n4), 1)); b4 = 0; end
                    n5 = find(VER_B.ID == TRI_B.vertices(T,5)); x5 = VER_B.x(n5); y5 = VER_B.y(n5); b5 = 1; if isempty(find(VER_B_boundary == VER_B.ID(n5), 1)); b5 = 0; end
                    n6 = find(VER_B.ID == TRI_B.vertices(T,6)); x6 = VER_B.x(n6); y6 = VER_B.y(n6); b6 = 1; if isempty(find(VER_B_boundary == VER_B.ID(n6), 1)); b6 = 0; end
                    if b4 == 1; WH(end+1) = 12; end
                    if b5 == 1; WH(end+1) = 23; end
                    if b6 == 1; WH(end+1) = 13; end
                    for W = 1:length(WH)
                        if WH(W) == 12; x_e = x4; y_e = y4; n = n4; N = 4; x_o1 = x5; y_o1 = y5; n_o1 = n5; N_o1 = 5; x_o2 = x6; y_o2 = y6; n_o2 = n6; N_o2 = 6; end
                        if WH(W) == 13; x_e = x6; y_e = y6; n = n6; N = 6; x_o1 = x4; y_o1 = y4; n_o1 = n4; N_o1 = 4; x_o2 = x5; y_o2 = y5; n_o2 = n5; N_o2 = 5; end
                        if WH(W) == 23; x_e = x5; y_e = y5; n = n5; N = 5; x_o1 = x4; y_o1 = y4; n_o1 = n4; N_o1 = 4; x_o2 = x6; y_o2 = y6; n_o2 = n6; N_o2 = 6; end
                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Ints2.x; OP.P.Ints2.y]);
                        PTx = OP.P.Ints2.x(IDX); PTy = OP.P.Ints2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Ints2.ID(IDX); ns(end+1) = n; end

                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Bots2.x; OP.P.Bots2.y]);
                        PTx = OP.P.Bots2.x(IDX); PTy = OP.P.Bots2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Bots2.ID(IDX); ns(end+1) = n; end

                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Lefts2.x; OP.P.Lefts2.y]);
                        PTx = OP.P.Lefts2.x(IDX); PTy = OP.P.Lefts2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Lefts2.ID(IDX); ns(end+1) = n; end

                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Rights2.x; OP.P.Rights2.y]);
                        PTx = OP.P.Rights2.x(IDX); PTy = OP.P.Rights2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Rights2.ID(IDX); ns(end+1) = n; end
                    end
                end
            end
            for T = 1:length(TRI_B.ID)
                for N = 1:length(OP.P.Bots.ID)
                    a = find(TRI_B.vertices(T,:) == OP.P.Bots.ID(N));
                    if ~ isempty(a); TRI_B.vertices(T,a) = OP.P.Tops.ID(N); end
                end
                for N2 = 1:length(OP.P.Bots2.ID)
                    a2 = find(TRI_B.vertices(T,:) == OP.P.Bots2.ID(N2));
                    if ~ isempty(a2); TRI_B.vertices(T,a2) = OP.P.Tops2.ID(N2); end
                end
            end
            if GEO.helix.order == 1; en = 1; end
            if GEO.helix.order == 2; VER_B.ID(ns) = []; VER_B.x(ns) = []; VER_B.y(ns) = [];
                if GEO.helix.element == 'j'; en = 3; else en = 2; end
            end

            ELEMENTS.nodes(end+1:end+size(TRI_B.vertices,1),:) = [TRI_B.vertices zeros(size(TRI_B.vertices,1),en)];
            ELEMENTS.index(end+1:end+length(TRI_B.ID)) = TRI_B.ID;
            if GEO.helix.order == 1; ELEMENTS.type(end+1:end+length(TRI_B.ID)) = char(ones(1,length(TRI_B.ID))*'d'); len = length(NODES.index); end
            if GEO.helix.order == 2; ELEMENTS.type(end+1:end+length(TRI_B.ID)) = char(ones(1,length(TRI_B.ID))*'f'); len = length(NODES.index); end
            NODES.index(end+1:end+length(VER_B.ID)) = VER_B.ID;
            NODES.x(end+1:end+length(VER_B.ID)) = VER_B.x;
            NODES.y(end+1:end+length(VER_B.ID)) = VER_B.y;
            for N = 1:length(VER_B.ID)
                % Radial coordinate
                NR = N_R(OP.R,d0,abs(NODES.y(len+N)-0.5*PPP),lambda);
                [NODES.X(len+N), NODES.Y(len+N), NODES.Z(len+N)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODES.x(len+N),NODES.y(len+N),NR);
            end
            OP.Hel_Tri = length(ELEMENTS.index) - OP.Hel_Quad;
        end

    case 2
        if isempty(find(master == rank,1)) == 0
            % CPU = 2 - Master Thread - Triangulating the LHS strip
            Nmax = max(NODES.index);
            MPI_Send(slave(1),1000000,PP.comm,Nmax,max(ELEMENTS.index),OP,GEO); % Sending info to slave thread about the RHS strip

            BL_x = []; BL_y = []; BL_ID = [];
            BL_x(1) = OP.P.LTop.x; BL_y(1) = OP.P.LTop.y; BL_ID(1) = max(NODES.index)+1;
            N = round(2*pi*OP.R/OP.D1);
            for Y = 1:N-1
                BL_x(end+1) = OP.P.LTop.x + Y*(OP.P.LBot.x - OP.P.LTop.x)/N;
                BL_y(end+1) = OP.P.LTop.y + Y*(OP.P.LBot.y - OP.P.LTop.y)/N;
                if OP.Tols.interior == 1
                    OP.P.Ls.x(end+1) = fL_x(BL_y(end),PPP,m) + (1/(OP.Tols.sideE+1))*OP.xoff;
                    OP.P.Ls.y(end+1) = BL_y(end);
                end
            end
            BL_ID(1:length(BL_x)) = max(BL_ID) + [1:length(BL_x)];
            BL_x(end+1) = OP.P.LBot.x; BL_y(end+1) = OP.P.LBot.y; BL_ID(end+1) = max(BL_ID) + 1;
            BL_x(end+1) = 0.5*(min(OP.P.Bots.x) + OP.P.LBot.x); BL_y(end+1) = OP.P.LBot.y; BL_ID(end+1) = max(BL_ID) + 1; bts(1:3) = [BL_ID(end-1:end) OP.P.Bots.ID(1)];
            BL_x(end+1:end+length(OP.P.Lefts.ID)) = OP.P.Lefts.x; BL_y(end+1:end+length(OP.P.Lefts.ID)) = OP.P.Lefts.y; BL_ID(end+1:end+length(OP.P.Lefts.ID)) = OP.P.Lefts.ID;
            BL_x(end+1) = 0.5*(min(OP.P.Tops.x) + OP.P.LTop.x); BL_y(end+1) = OP.P.LTop.y; BL_ID(end+1) = max(BL_ID) + 1; tps(1:3) = [BL_ID(1) BL_ID(end) OP.P.Tops.ID(1)];
            [BL_ID, ids] = EXTERNAL_unique_no_sort(BL_ID); BL_x = BL_x(ids); BL_y = BL_y(ids);
            if OP.Tols.interior == 1
                IL_x = OP.P.Ls.x; IL_y = OP.P.Ls.y; IL_ID = [];
                IL_ID(1:length(IL_x)) = max(BL_ID) + [1:length(IL_x)];
            else IL_x = []; IL_y = []; IL_ID = []; end
            
            [TRI_L,VER_L,TRI_L_boundary,VER_L_boundary] = AJSABQ_Delaunay(BL_x,BL_y,BL_ID,IL_x,IL_y,IL_ID,OP.Tols.tolRND,OP.Tols.tolNN,OP.Tols.tolBG,GEO.helix.order,max(NODES.index)+1,max(ELEMENTS.index)+1,OP.Tols.scale);
            disp('Consolidating');
            if GEO.helix.order == 2
                ns = []; bt = find(VER_L.y < min(VER_L.y) + OP.Tols.tolNN); tp = find(VER_L.y > max(VER_L.y) - OP.Tols.tolNN); 
                bts(end+1:end+length(bt)) = VER_L.ID(bt); tps(end+1:end+length(tp)) = VER_L.ID(tp);
                bts = unique(bts); tps = unique(tps); bts_x = []; tps_x = [];
                for bt = 1:length(bts)
                    bts_x(bt) = VER_L.x(find(VER_L.ID == bts(bt)));
                    tps_x(bt) = VER_L.x(find(VER_L.ID == tps(bt)));
                end
                [bts_x, b] = sort(bts_x); [tps_x, t] = sort(tps_x);
                bts = bts(b); tps = tps(t);
            end
            for T = 1:length(TRI_L.ID)
                n1 = find(VER_L.ID == TRI_L.vertices(T,1)); x1 = VER_L.x(n1); y1 = VER_L.y(n1);
                n2 = find(VER_L.ID == TRI_L.vertices(T,2)); x2 = VER_L.x(n2); y2 = VER_L.y(n2);
                n3 = find(VER_L.ID == TRI_L.vertices(T,3)); x3 = VER_L.x(n3); y3 = VER_L.y(n3);
                if GEO.helix.order == 2; WH = [];
                    n4 = find(VER_L.ID == TRI_L.vertices(T,4)); x4 = VER_L.x(n4); y4 = VER_L.y(n4); b4 = 1; if isempty(find(VER_L_boundary == VER_L.ID(n4), 1)); b4 = 0; end
                    n5 = find(VER_L.ID == TRI_L.vertices(T,5)); x5 = VER_L.x(n5); y5 = VER_L.y(n5); b5 = 1; if isempty(find(VER_L_boundary == VER_L.ID(n5), 1)); b5 = 0; end
                    n6 = find(VER_L.ID == TRI_L.vertices(T,6)); x6 = VER_L.x(n6); y6 = VER_L.y(n6); b6 = 1; if isempty(find(VER_L_boundary == VER_L.ID(n6), 1)); b6 = 0; end
                    if b4 == 1; WH(end+1) = 12; end
                    if b5 == 1; WH(end+1) = 23; end
                    if b6 == 1; WH(end+1) = 13; end
                    for W = 1:length(WH)
                        if WH(W) == 12; x_e = x4; y_e = y4; n = n4; N = 4; end
                        if WH(W) == 13; x_e = x6; y_e = y6; n = n6; N = 6; end
                        if WH(W) == 23; x_e = x5; y_e = y5; n = n5; N = 5; end
                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Lefts2.x; OP.P.Lefts2.y]);
                        PTx = OP.P.Lefts2.x(IDX); PTy = OP.P.Lefts2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_L.vertices(T,N) = OP.P.Lefts2.ID(IDX); ns(end+1) = n; end
                    end
                end
                x_edge12 = 0.5*(x1 + x2); y_edge12 = 0.5*(y1 + y2);
                x_edge13 = 0.5*(x1 + x3); y_edge13 = 0.5*(y1 + y3);
                x_edge23 = 0.5*(x2 + x3); y_edge23 = 0.5*(y2 + y3);
                if abs(x_edge12 - fL_x(y_edge12,PPP,m)) < OP.Tols.tolNN
                    NSET.Left(end+1) = VER_L.ID(n1); NSET.Left(end+1) = VER_L.ID(n2); ELSET.Left(end+1) = TRI_L.ID(T);
                    if GEO.helix.order == 1; TRI_L.vertices(T,:) = [VER_L.ID(n2) VER_L.ID(n3) VER_L.ID(n1)]; end
                    if GEO.helix.order == 2; TRI_L.vertices(T,:) = [VER_L.ID(n2) VER_L.ID(n3) VER_L.ID(n1) VER_L.ID(n5) VER_L.ID(n6) VER_L.ID(n4)]; NSET.Left(end+1) = VER_L.ID(n4); end
                end
                if abs(x_edge13 - fL_x(y_edge13,PPP,m)) < OP.Tols.tolNN
                    NSET.Left(end+1) = VER_L.ID(n1); NSET.Left(end+1) = VER_L.ID(n3); ELSET.Left(end+1) = TRI_L.ID(T);
                    if GEO.helix.order == 2; NSET.Left(end+1) = VER_L.ID(n6); end
                end
                if abs(x_edge23 - fL_x(y_edge23,PPP,m)) < OP.Tols.tolNN
                    NSET.Left(end+1) = VER_L.ID(n2); NSET.Left(end+1) = VER_L.ID(n3); ELSET.Left(end+1) = TRI_L.ID(T);
                    if GEO.helix.order == 1; TRI_L.vertices(T,:) = [VER_L.ID(n3) VER_L.ID(n1) VER_L.ID(n2)]; end
                    if GEO.helix.order == 2; TRI_L.vertices(T,:) = [VER_L.ID(n3) VER_L.ID(n1) VER_L.ID(n2) VER_L.ID(n6) VER_L.ID(n4) VER_L.ID(n5)]; NSET.Left(end+1) = VER_L.ID(n5); end
                end
                for L = 1:length(bts)
                    n = find(TRI_L.vertices(T,:) == bts(L));
                    if ~ isempty(n); TRI_L.vertices(T,n) = tps(L); end
                end
                if (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) < 0.0 % Verifying the consistency of the element normal
                   if GEO.helix.order == 1; TRI_L.vertices(T,:) = [TRI_L.vertices(T,3) TRI_L.vertices(T,2) TRI_L.vertices(T,1)]; end
                   if GEO.helix.order == 2; TRI_L.vertices(T,:) = [TRI_L.vertices(T,3) TRI_L.vertices(T,2) TRI_L.vertices(T,1) TRI_L.vertices(T,5) TRI_L.vertices(T,4) TRI_L.vertices(T,6)]; end                   
                end
            end
            if GEO.helix.order == 1; en = 1; end
            if GEO.helix.order == 2; VER_L.ID(ns) = []; VER_L.x(ns) = []; VER_L.y(ns) = [];
                if GEO.helix.element == 'j'; en = 3; else en = 2; end
            end
            NSET.Bottom = unique(NSET.Left); NSET = rmfield(NSET,'Left'); NSET = rmfield(NSET,'Right');
            ELSET.Bottom = ELSET.Left; ELSET = rmfield(ELSET,'Left'); ELSET = rmfield(ELSET,'Right');

            ELEMENTS.nodes(end+1:end+size(TRI_L.vertices,1),:) = [TRI_L.vertices zeros(size(TRI_L.vertices,1),en)];
            ELEMENTS.index(end+1:end+length(TRI_L.ID)) = TRI_L.ID;
            if GEO.helix.order == 1; ELEMENTS.type(end+1:end+length(TRI_L.ID)) = char(ones(1,length(TRI_L.ID))*'d'); len = length(NODES.index); end
            if GEO.helix.order == 2; ELEMENTS.type(end+1:end+length(TRI_L.ID)) = char(ones(1,length(TRI_L.ID))*'f'); len = length(NODES.index); end
            NODES.index(end+1:end+length(VER_L.ID)) = VER_L.ID;
            NODES.x(end+1:end+length(VER_L.ID)) = VER_L.x;
            NODES.y(end+1:end+length(VER_L.ID)) = VER_L.y;
            for N = 1:length(VER_L.ID)
                % Radial coordinate
                NR = N_R(OP.R,d0,abs(NODES.y(len+N)-0.5*PPP),lambda);
                [NODES.X(len+N), NODES.Y(len+N), NODES.Z(len+N)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODES.x(len+N),NODES.y(len+N),NR);
            end


            % CPU = 2 - Master Thread - Triangulating the misalingment
            BB_ID = [OP.P.Bots.ID OP.P.Ints.ID]; BB_x = [OP.P.Bots.x OP.P.Ints.x]; BB_y = [OP.P.Bots.y OP.P.Ints.y]; ns = [];

            [TRI_B,VER_B,TRI_B_boundary,VER_B_boundary] = AJSABQ_Delaunay(BB_x,BB_y,BB_ID,[],[],[],OP.Tols.tolRND,OP.Tols.tolNN,OP.Tols.tolBG,GEO.helix.order,max(NODES.index)+1,max(ELEMENTS.index)+1,0);
            disp('Consolidating');
            for T = 1:length(TRI_B.ID)
                if GEO.helix.order == 2; WH = [];
                    n4 = find(VER_B.ID == TRI_B.vertices(T,4)); x4 = VER_B.x(n4); y4 = VER_B.y(n4); b4 = 1; if isempty(find(VER_B_boundary == VER_B.ID(n4), 1)); b4 = 0; end
                    n5 = find(VER_B.ID == TRI_B.vertices(T,5)); x5 = VER_B.x(n5); y5 = VER_B.y(n5); b5 = 1; if isempty(find(VER_B_boundary == VER_B.ID(n5), 1)); b5 = 0; end
                    n6 = find(VER_B.ID == TRI_B.vertices(T,6)); x6 = VER_B.x(n6); y6 = VER_B.y(n6); b6 = 1; if isempty(find(VER_B_boundary == VER_B.ID(n6), 1)); b6 = 0; end
                    if b4 == 1; WH(end+1) = 12; end
                    if b5 == 1; WH(end+1) = 23; end
                    if b6 == 1; WH(end+1) = 13; end
                    for W = 1:length(WH)
                        if WH(W) == 12; x_e = x4; y_e = y4; n = n4; N = 4; x_o1 = x5; y_o1 = y5; n_o1 = n5; N_o1 = 5; x_o2 = x6; y_o2 = y6; n_o2 = n6; N_o2 = 6; end
                        if WH(W) == 13; x_e = x6; y_e = y6; n = n6; N = 6; x_o1 = x4; y_o1 = y4; n_o1 = n4; N_o1 = 4; x_o2 = x5; y_o2 = y5; n_o2 = n5; N_o2 = 5; end
                        if WH(W) == 23; x_e = x5; y_e = y5; n = n5; N = 5; x_o1 = x4; y_o1 = y4; n_o1 = n4; N_o1 = 4; x_o2 = x6; y_o2 = y6; n_o2 = n6; N_o2 = 6; end
                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Ints2.x; OP.P.Ints2.y]);
                        PTx = OP.P.Ints2.x(IDX); PTy = OP.P.Ints2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Ints2.ID(IDX); ns(end+1) = n; end

                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Bots2.x; OP.P.Bots2.y]);
                        PTx = OP.P.Bots2.x(IDX); PTy = OP.P.Bots2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Bots2.ID(IDX); ns(end+1) = n; end

                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Lefts2.x; OP.P.Lefts2.y]);
                        PTx = OP.P.Lefts2.x(IDX); PTy = OP.P.Lefts2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Lefts2.ID(IDX); ns(end+1) = n; end

                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Rights2.x; OP.P.Rights2.y]);
                        PTx = OP.P.Rights2.x(IDX); PTy = OP.P.Rights2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Rights2.ID(IDX); ns(end+1) = n; end
                    end
                end
            end
            for T = 1:length(TRI_B.ID)
                for N = 1:length(OP.P.Bots.ID)
                    a = find(TRI_B.vertices(T,:) == OP.P.Bots.ID(N));
                    if ~ isempty(a); TRI_B.vertices(T,a) = OP.P.Tops.ID(N); end
                end
                for N2 = 1:length(OP.P.Bots2.ID)
                    a2 = find(TRI_B.vertices(T,:) == OP.P.Bots2.ID(N2));
                    if ~ isempty(a2); TRI_B.vertices(T,a2) = OP.P.Tops2.ID(N2); end
                end
            end
            if GEO.helix.order == 1; en = 1; end
            if GEO.helix.order == 2; VER_B.ID(ns) = []; VER_B.x(ns) = []; VER_B.y(ns) = [];
                if GEO.helix.element == 'j'; en = 3; else en = 2; end
            end

            ELEMENTS.nodes(end+1:end+size(TRI_B.vertices,1),:) = [TRI_B.vertices zeros(size(TRI_B.vertices,1),en)];
            ELEMENTS.index(end+1:end+length(TRI_B.ID)) = TRI_B.ID;
            if GEO.helix.order == 1; ELEMENTS.type(end+1:end+length(TRI_B.ID)) = char(ones(1,length(TRI_B.ID))*'d'); len = length(NODES.index); end
            if GEO.helix.order == 2; ELEMENTS.type(end+1:end+length(TRI_B.ID)) = char(ones(1,length(TRI_B.ID))*'f'); len = length(NODES.index); end
            NODES.index(end+1:end+length(VER_B.ID)) = VER_B.ID;
            NODES.x(end+1:end+length(VER_B.ID)) = VER_B.x;
            NODES.y(end+1:end+length(VER_B.ID)) = VER_B.y;
            for N = 1:length(VER_B.ID)
                % Radial coordinate
                NR = N_R(OP.R,d0,abs(NODES.y(len+N)-0.5*PPP),lambda);
                [NODES.X(len+N), NODES.Y(len+N), NODES.Z(len+N)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODES.x(len+N),NODES.y(len+N),NR);
            end

        elseif isempty(find(slave == rank,1)) == 0 && PP.rank == 1
            [Nmax, Emax, OP, GEO] = MPI_Recv(master,1000000,PP.comm); % Receiving info from master thread about the RHS strip
            NODES = struct('index',[],'x',[],'y',[],'X',[],'Y',[],'Z',[],'case',[]);
            ELEMENTS = struct('index',[],'nodes',[],'type',[],'case',[]);
            NSET = struct('Top',[],'Right',[]);
            ELSET = struct('Top',[],'Right',[]);

            % CPU = 2 - Slave 1 Thread - Triangulating the RHS strip
            BR_x = []; BR_y = []; BR_ID = []; clear('bts','tps');
            BR_x(1) = OP.P.RTop.x; BR_y(1) = OP.P.RTop.y; BR_ID(1) = Nmax*2+1;
            N = round(2*pi*OP.R/OP.D1);
            for Y = 1:N-1
                BR_x(end+1) = OP.P.RTop.x + Y*(OP.P.RBot.x - OP.P.RTop.x)/N;
                BR_y(end+1) = OP.P.RTop.y + Y*(OP.P.RBot.y - OP.P.RTop.y)/N;
                if OP.Tols.interior == 1
                    OP.P.Rs.x(end+1) = fR_x(BR_y(end),x_max,PPP,m) - (1/(OP.Tols.sideE+1))*OP.xoff;
                    OP.P.Rs.y(end+1) = BR_y(end);
                end
            end
            BR_ID(1:length(BR_x)) = max(BR_ID) + [1:length(BR_x)];
            BR_x(end+1) = OP.P.RBot.x; BR_y(end+1) = OP.P.RBot.y; BR_ID(end+1) = max(BR_ID) + 1;
            BR_x(end+1) = 0.5*(max(OP.P.Bots.x) + OP.P.RBot.x); BR_y(end+1) = OP.P.RBot.y; BR_ID(end+1) = max(BR_ID) + 1; bts(1:3) = [BR_ID(end-1:end) OP.P.Bots.ID(end)]; bts = fliplr(bts);
            BR_x(end+1:end+length(OP.P.Rights.ID)) = OP.P.Rights.x; BR_y(end+1:end+length(OP.P.Rights.ID)) = OP.P.Rights.y; BR_ID(end+1:end+length(OP.P.Rights.ID)) = OP.P.Rights.ID;
            BR_x(end+1) = 0.5*(max(OP.P.Tops.x) + OP.P.RTop.x); BR_y(end+1) = OP.P.RTop.y; BR_ID(end+1) = max(BR_ID) + 1; tps(1:3) = [BR_ID(1) BR_ID(end) OP.P.Tops.ID(end)]; tps = fliplr(tps);
            [BR_ID, ids] = EXTERNAL_unique_no_sort(BR_ID); BR_x = BR_x(ids); BR_y = BR_y(ids);
            if OP.Tols.interior == 1
                IR_x = OP.P.Rs.x; IR_y = OP.P.Rs.y; IR_ID = [];
                IR_ID(1:length(IR_x)) = max(BR_ID) + [1:length(IR_x)];
            else IR_x = []; IR_y = []; IR_ID = []; end
            
            [TRI_R,VER_R,TRI_R_boundary,VER_R_boundary] = AJSABQ_Delaunay(BR_x,BR_y,BR_ID,IR_x,IR_y,IR_ID,OP.Tols.tolRND,OP.Tols.tolNN,OP.Tols.tolBG,GEO.helix.order,Nmax*2+1,Emax*2+1,OP.Tols.scale);
            disp('Consolidating');
            if GEO.helix.order == 2
                ns = []; bt = find(VER_R.y < min(VER_R.y) + OP.Tols.tolNN); tp = find(VER_R.y > max(VER_R.y) - OP.Tols.tolNN); 
                bts(end+1:end+length(bt)) = VER_R.ID(bt); tps(end+1:end+length(tp)) = VER_R.ID(tp);
                bts = unique(bts); tps = unique(tps); bts_x = []; tps_x = [];
                for bt = 1:length(bts)
                    bts_x(bt) = VER_R.x(find(VER_R.ID == bts(bt)));
                    tps_x(bt) = VER_R.x(find(VER_R.ID == tps(bt)));
                end
                [bts_x, b] = sort(bts_x); [tps_x, t] = sort(tps_x);
                bts = bts(b); tps = tps(t);
            end
            for T = 1:length(TRI_R.ID)
                n1 = find(VER_R.ID == TRI_R.vertices(T,1)); x1 = VER_R.x(n1); y1 = VER_R.y(n1);
                n2 = find(VER_R.ID == TRI_R.vertices(T,2)); x2 = VER_R.x(n2); y2 = VER_R.y(n2);
                n3 = find(VER_R.ID == TRI_R.vertices(T,3)); x3 = VER_R.x(n3); y3 = VER_R.y(n3);
                if GEO.helix.order == 2; WH = [];
                    n4 = find(VER_R.ID == TRI_R.vertices(T,4)); x4 = VER_R.x(n4); y4 = VER_R.y(n4); b4 = 1; if isempty(find(VER_R_boundary == VER_R.ID(n4), 1)); b4 = 0; end
                    n5 = find(VER_R.ID == TRI_R.vertices(T,5)); x5 = VER_R.x(n5); y5 = VER_R.y(n5); b5 = 1; if isempty(find(VER_R_boundary == VER_R.ID(n5), 1)); b5 = 0; end
                    n6 = find(VER_R.ID == TRI_R.vertices(T,6)); x6 = VER_R.x(n6); y6 = VER_R.y(n6); b6 = 1; if isempty(find(VER_R_boundary == VER_R.ID(n6), 1)); b6 = 0; end
                    if b4 == 1; WH(end+1) = 12; end
                    if b5 == 1; WH(end+1) = 23; end
                    if b6 == 1; WH(end+1) = 13; end
                    for W = 1:length(WH)
                        if WH(W) == 12; x_e = x4; y_e = y4; n = n4; N = 4; end
                        if WH(W) == 13; x_e = x6; y_e = y6; n = n6; N = 6; end
                        if WH(W) == 23; x_e = x5; y_e = y5; n = n5; N = 5; end
                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Rights2.x; OP.P.Rights2.y]);
                        PTx = OP.P.Rights2.x(IDX); PTy = OP.P.Rights2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_R.vertices(T,N) = OP.P.Rights2.ID(IDX); ns(end+1) = n; end
                    end
                end
                x_edge12 = 0.5*(x1 + x2); y_edge12 = 0.5*(y1 + y2);
                x_edge13 = 0.5*(x1 + x3); y_edge13 = 0.5*(y1 + y3);
                x_edge23 = 0.5*(x2 + x3); y_edge23 = 0.5*(y2 + y3);
                if abs(x_edge12 - fR_x(y_edge12,x_max,PPP,m)) < OP.Tols.tolNN
                    NSET.Right(end+1) = VER_R.ID(n1); NSET.Right(end+1) = VER_R.ID(n2); ELSET.Right(end+1) = TRI_R.ID(T);
                    if GEO.helix.order == 1; TRI_R.vertices(T,:) = [VER_R.ID(n2) VER_R.ID(n3) VER_R.ID(n1)]; end
                    if GEO.helix.order == 2; TRI_R.vertices(T,:) = [VER_R.ID(n2) VER_R.ID(n3) VER_R.ID(n1) VER_R.ID(n5) VER_R.ID(n6) VER_R.ID(n4)]; NSET.Right(end+1) = VER_R.ID(n4); end
                end
                if abs(x_edge13 - fR_x(y_edge13,x_max,PPP,m)) < OP.Tols.tolNN
                    NSET.Right(end+1) = VER_R.ID(n1); NSET.Right(end+1) = VER_R.ID(n3); ELSET.Right(end+1) = TRI_R.ID(T);
                    if GEO.helix.order == 2; NSET.Right(end+1) = VER_R.ID(n6); end
                end
                if abs(x_edge23 - fR_x(y_edge23,x_max,PPP,m)) < OP.Tols.tolNN
                    NSET.Right(end+1) = VER_R.ID(n2); NSET.Right(end+1) = VER_R.ID(n3); ELSET.Right(end+1) = TRI_R.ID(T);
                    if GEO.helix.order == 1; TRI_R.vertices(T,:) = [VER_R.ID(n3) VER_R.ID(n1) VER_R.ID(n2)]; end
                    if GEO.helix.order == 2; TRI_R.vertices(T,:) = [VER_R.ID(n3) VER_R.ID(n1) VER_R.ID(n2) VER_R.ID(n6) VER_R.ID(n4) VER_R.ID(n5)]; NSET.Right(end+1) = VER_R.ID(n5); end
                end
                for R = 1:length(bts)
                    n = find(TRI_R.vertices(T,:) == bts(R));
                    if ~ isempty(n); TRI_R.vertices(T,n) = tps(R); end
                end
                if (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) < 0.0 % Verifying the consistency of the element normal
                   if GEO.helix.order == 1; TRI_R.vertices(T,:) = [TRI_R.vertices(T,3) TRI_R.vertices(T,2) TRI_R.vertices(T,1)]; end
                   if GEO.helix.order == 2; TRI_R.vertices(T,:) = [TRI_R.vertices(T,3) TRI_R.vertices(T,2) TRI_R.vertices(T,1) TRI_R.vertices(T,5) TRI_R.vertices(T,4) TRI_R.vertices(T,6)]; end                   
                end                
            end
            if GEO.helix.order == 1; en = 1; end
            if GEO.helix.order == 2
                VER_R.ID(ns) = []; VER_R.x(ns) = []; VER_R.y(ns) = [];
                if GEO.helix.element == 'j'; en = 3; else en = 2; end
            end
            NSET.Top = unique(NSET.Right); NSET = rmfield(NSET,'Right');
            ELSET.Top = ELSET.Right; ELSET = rmfield(ELSET,'Right');

            ELEMENTS.nodes(end+1:end+size(TRI_R.vertices,1),:) = [TRI_R.vertices zeros(size(TRI_R.vertices,1),en)];
            ELEMENTS.index(end+1:end+length(TRI_R.ID)) = TRI_R.ID;
            if GEO.helix.order == 1; ELEMENTS.type(end+1:end+length(TRI_R.ID)) = char(ones(1,length(TRI_R.ID))*'d'); len = length(NODES.index); end
            if GEO.helix.order == 2; ELEMENTS.type(end+1:end+length(TRI_R.ID)) = char(ones(1,length(TRI_R.ID))*'f'); len = length(NODES.index); end
            NODES.index(end+1:end+length(VER_R.ID)) = VER_R.ID;
            NODES.x(end+1:end+length(VER_R.ID)) = VER_R.x;
            NODES.y(end+1:end+length(VER_R.ID)) = VER_R.y;
            for N = 1:length(VER_R.ID)
                % Radial coordinate
                NR = N_R(OP.R,d0,abs(NODES.y(len+N)-0.5*PPP),lambda);
                [NODES.X(len+N), NODES.Y(len+N), NODES.Z(len+N)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODES.x(len+N),NODES.y(len+N),NR);
            end

            MPI_Send(master,2000000,PP.comm,ELEMENTS,NODES,NSET.Top,ELSET.Top);
        end

        % Calibrating
        if isempty(find(master == rank,1)) == 0
            [SlaveELEMENTS, SlaveNODES, SlaveNTop, SlaveETop] = MPI_Recv(slave(1),2000000,PP.comm);
            ELEMENTS.nodes(end+1:end+length(SlaveELEMENTS.index),:) = SlaveELEMENTS.nodes;
            ELEMENTS.index(end+1:end+length(SlaveELEMENTS.index)) = SlaveELEMENTS.index;
            NODES.index(end+1:end+length(SlaveNODES.index)) = SlaveNODES.index;
            NODES.x(end+1:end+length(SlaveNODES.x)) = SlaveNODES.x;
            NODES.y(end+1:end+length(SlaveNODES.y)) = SlaveNODES.y;
            NODES.X(end+1:end+length(SlaveNODES.X)) = SlaveNODES.X;
            NODES.Y(end+1:end+length(SlaveNODES.Y)) = SlaveNODES.Y;
            NODES.Z(end+1:end+length(SlaveNODES.Z)) = SlaveNODES.Z;
            NSET.Top = SlaveNTop; ELSET.Top = SlaveETop;
            OP.Hel_Tri = length(ELEMENTS.index) - OP.Hel_Quad;
        end

    otherwise
        if isempty(find(master == rank,1)) == 0           
            % CPU >= 3 - Master Thread - Triangulating the LHS strip
            Nmax = max(NODES.index);
            MPI_Send(slave(1),1000000,PP.comm,Nmax,max(ELEMENTS.index),OP,GEO); % Sending info to slave thread about the RHS strip
            MPI_Send(slave(2),1500000,PP.comm,Nmax,max(ELEMENTS.index),OP,GEO); % Sending info to slave thread about the misalingment

            BL_x = []; BL_y = []; BL_ID = [];
            BL_x(1) = OP.P.LTop.x; BL_y(1) = OP.P.LTop.y; BL_ID(1) = max(NODES.index)+1;
            N = round(2*pi*OP.R/OP.D1);
            for Y = 1:N-1
                BL_x(end+1) = OP.P.LTop.x + Y*(OP.P.LBot.x - OP.P.LTop.x)/N;
                BL_y(end+1) = OP.P.LTop.y + Y*(OP.P.LBot.y - OP.P.LTop.y)/N;
                if OP.Tols.interior == 1
                    OP.P.Ls.x(end+1) = fL_x(BL_y(end),PPP,m) + (1/(OP.Tols.sideE+1))*OP.xoff;
                    OP.P.Ls.y(end+1) = BL_y(end);
                end
            end
            BL_ID(1:length(BL_x)) = max(BL_ID) + [1:length(BL_x)];
            BL_x(end+1) = OP.P.LBot.x; BL_y(end+1) = OP.P.LBot.y; BL_ID(end+1) = max(BL_ID) + 1;
            BL_x(end+1) = 0.5*(min(OP.P.Bots.x) + OP.P.LBot.x); BL_y(end+1) = OP.P.LBot.y; BL_ID(end+1) = max(BL_ID) + 1; bts(1:3) = [BL_ID(end-1:end) OP.P.Bots.ID(1)];
            BL_x(end+1:end+length(OP.P.Lefts.ID)) = OP.P.Lefts.x; BL_y(end+1:end+length(OP.P.Lefts.ID)) = OP.P.Lefts.y; BL_ID(end+1:end+length(OP.P.Lefts.ID)) = OP.P.Lefts.ID;
            BL_x(end+1) = 0.5*(min(OP.P.Tops.x) + OP.P.LTop.x); BL_y(end+1) = OP.P.LTop.y; BL_ID(end+1) = max(BL_ID) + 1; tps(1:3) = [BL_ID(1) BL_ID(end) OP.P.Tops.ID(1)];
            [BL_ID, ids] = EXTERNAL_unique_no_sort(BL_ID); BL_x = BL_x(ids); BL_y = BL_y(ids);
            if OP.Tols.interior == 1
                IL_x = OP.P.Ls.x; IL_y = OP.P.Ls.y; IL_ID = [];
                IL_ID(1:length(IL_x)) = max(BL_ID) + [1:length(IL_x)];
            else IL_x = []; IL_y = []; IL_ID = []; end
            
            [TRI_L,VER_L,TRI_L_boundary,VER_L_boundary] = AJSABQ_Delaunay(BL_x,BL_y,BL_ID,IL_x,IL_y,IL_ID,OP.Tols.tolRND,OP.Tols.tolNN,OP.Tols.tolBG,GEO.helix.order,max(NODES.index)+1,max(ELEMENTS.index)+1,OP.Tols.scale);
            disp('Consolidating');
            if GEO.helix.order == 2
                ns = []; bt = find(VER_L.y < min(VER_L.y) + OP.Tols.tolNN); tp = find(VER_L.y > max(VER_L.y) - OP.Tols.tolNN); 
                bts(end+1:end+length(bt)) = VER_L.ID(bt); tps(end+1:end+length(tp)) = VER_L.ID(tp);
                bts = unique(bts); tps = unique(tps); bts_x = []; tps_x = [];  
                for bt = 1:length(bts);
                    bts_x(bt) = VER_L.x(find(VER_L.ID == bts(bt)));
                    tps_x(bt) = VER_L.x(find(VER_L.ID == tps(bt)));
                end
                [bts_x, b] = sort(bts_x); [tps_x, t] = sort(tps_x);
                bts = bts(b); tps = tps(t);
            end
            for T = 1:length(TRI_L.ID)
                n1 = find(VER_L.ID == TRI_L.vertices(T,1)); x1 = VER_L.x(n1); y1 = VER_L.y(n1);
                n2 = find(VER_L.ID == TRI_L.vertices(T,2)); x2 = VER_L.x(n2); y2 = VER_L.y(n2);
                n3 = find(VER_L.ID == TRI_L.vertices(T,3)); x3 = VER_L.x(n3); y3 = VER_L.y(n3);
                if GEO.helix.order == 2; WH = [];
                    n4 = find(VER_L.ID == TRI_L.vertices(T,4)); x4 = VER_L.x(n4); y4 = VER_L.y(n4); b4 = 1; if isempty(find(VER_L_boundary == VER_L.ID(n4), 1)); b4 = 0; end
                    n5 = find(VER_L.ID == TRI_L.vertices(T,5)); x5 = VER_L.x(n5); y5 = VER_L.y(n5); b5 = 1; if isempty(find(VER_L_boundary == VER_L.ID(n5), 1)); b5 = 0; end
                    n6 = find(VER_L.ID == TRI_L.vertices(T,6)); x6 = VER_L.x(n6); y6 = VER_L.y(n6); b6 = 1; if isempty(find(VER_L_boundary == VER_L.ID(n6), 1)); b6 = 0; end
                    if b4 == 1; WH(end+1) = 12; end
                    if b5 == 1; WH(end+1) = 23; end
                    if b6 == 1; WH(end+1) = 13; end
                    for W = 1:length(WH)
                        if WH(W) == 12; x_e = x4; y_e = y4; n = n4; N = 4; end
                        if WH(W) == 13; x_e = x6; y_e = y6; n = n6; N = 6; end
                        if WH(W) == 23; x_e = x5; y_e = y5; n = n5; N = 5; end
                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Lefts2.x; OP.P.Lefts2.y]);
                        PTx = OP.P.Lefts2.x(IDX); PTy = OP.P.Lefts2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_L.vertices(T,N) = OP.P.Lefts2.ID(IDX); ns(end+1) = n; end
                    end
                end
                x_edge12 = 0.5*(x1 + x2); y_edge12 = 0.5*(y1 + y2);
                x_edge13 = 0.5*(x1 + x3); y_edge13 = 0.5*(y1 + y3);
                x_edge23 = 0.5*(x2 + x3); y_edge23 = 0.5*(y2 + y3);
                if abs(x_edge12 - fL_x(y_edge12,PPP,m)) < OP.Tols.tolNN
                    NSET.Left(end+1) = VER_L.ID(n1); NSET.Left(end+1) = VER_L.ID(n2); ELSET.Left(end+1) = TRI_L.ID(T);
                    if GEO.helix.order == 1; TRI_L.vertices(T,:) = [VER_L.ID(n2) VER_L.ID(n3) VER_L.ID(n1)]; end
                    if GEO.helix.order == 2; TRI_L.vertices(T,:) = [VER_L.ID(n2) VER_L.ID(n3) VER_L.ID(n1) VER_L.ID(n5) VER_L.ID(n6) VER_L.ID(n4)]; NSET.Left(end+1) = VER_L.ID(n4); end
                end
                if abs(x_edge13 - fL_x(y_edge13,PPP,m)) < OP.Tols.tolNN
                    NSET.Left(end+1) = VER_L.ID(n1); NSET.Left(end+1) = VER_L.ID(n3); ELSET.Left(end+1) = TRI_L.ID(T);
                    if GEO.helix.order == 2; NSET.Left(end+1) = VER_L.ID(n6); end
                end
                if abs(x_edge23 - fL_x(y_edge23,PPP,m)) < OP.Tols.tolNN
                    NSET.Left(end+1) = VER_L.ID(n2); NSET.Left(end+1) = VER_L.ID(n3); ELSET.Left(end+1) = TRI_L.ID(T);
                    if GEO.helix.order == 1; TRI_L.vertices(T,:) = [VER_L.ID(n3) VER_L.ID(n1) VER_L.ID(n2)]; end
                    if GEO.helix.order == 2; TRI_L.vertices(T,:) = [VER_L.ID(n3) VER_L.ID(n1) VER_L.ID(n2) VER_L.ID(n6) VER_L.ID(n4) VER_L.ID(n5)]; NSET.Left(end+1) = VER_L.ID(n5); end
                end
                for L = 1:length(bts)
                    n = find(TRI_L.vertices(T,:) == bts(L));
                    if ~ isempty(n); TRI_L.vertices(T,n) = tps(L); end
                end 
                if (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) < 0.0 % Verifying the consistency of the element normal
                   if GEO.helix.order == 1; TRI_L.vertices(T,:) = [TRI_L.vertices(T,3) TRI_L.vertices(T,2) TRI_L.vertices(T,1)]; end
                   if GEO.helix.order == 2; TRI_L.vertices(T,:) = [TRI_L.vertices(T,3) TRI_L.vertices(T,2) TRI_L.vertices(T,1) TRI_L.vertices(T,5) TRI_L.vertices(T,4) TRI_L.vertices(T,6)]; end                   
                end
            end
            if GEO.helix.order == 1; en = 1; end
            if GEO.helix.order == 2; VER_L.ID(ns) = []; VER_L.x(ns) = []; VER_L.y(ns) = [];
                if GEO.helix.element == 'j'; en = 3; else en = 2; end
            end
            NSET.Bottom = unique(NSET.Left); NSET = rmfield(NSET,'Left'); NSET = rmfield(NSET,'Right');
            ELSET.Bottom = ELSET.Left; ELSET = rmfield(ELSET,'Left'); ELSET = rmfield(ELSET,'Right');

            ELEMENTS.nodes(end+1:end+size(TRI_L.vertices,1),:) = [TRI_L.vertices zeros(size(TRI_L.vertices,1),en)];
            ELEMENTS.index(end+1:end+length(TRI_L.ID)) = TRI_L.ID;
            if GEO.helix.order == 1; ELEMENTS.type(end+1:end+length(TRI_L.ID)) = char(ones(1,length(TRI_L.ID))*'d'); len = length(NODES.index); end
            if GEO.helix.order == 2; ELEMENTS.type(end+1:end+length(TRI_L.ID)) = char(ones(1,length(TRI_L.ID))*'f'); len = length(NODES.index); end
            NODES.index(end+1:end+length(VER_L.ID)) = VER_L.ID;
            NODES.x(end+1:end+length(VER_L.ID)) = VER_L.x;
            NODES.y(end+1:end+length(VER_L.ID)) = VER_L.y;
            for N = 1:length(VER_L.ID)
                % Radial coordinate
                NR = N_R(OP.R,d0,abs(NODES.y(len+N)-0.5*PPP),lambda);
                [NODES.X(len+N), NODES.Y(len+N), NODES.Z(len+N)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODES.x(len+N),NODES.y(len+N),NR);
            end

        elseif isempty(find(slave == rank,1)) == 0 && rank == 1
            [Nmax, Emax, OP, GEO] = MPI_Recv(master,1000000,PP.comm); % Receiving info from master thread about the RHS strip

            NODES = struct('index',[],'x',[],'y',[],'X',[],'Y',[],'Z',[],'case',[]);
            ELEMENTS = struct('index',[],'nodes',[],'type',[],'case',[]);
            NSET = struct('Top',[],'Right',[]);
            ELSET = struct('Top',[],'Right',[]);

            % CPU >= 3 - Slave 1 Thread - Triangulating the RHS strip
            BR_x = []; BR_y = []; BR_ID = []; clear('bts','tps');
            BR_x(1) = OP.P.RTop.x; BR_y(1) = OP.P.RTop.y; BR_ID(1) = Nmax*2+1;
            N = round(2*pi*OP.R/OP.D1);
            for Y = 1:N-1
                BR_x(end+1) = OP.P.RTop.x + Y*(OP.P.RBot.x - OP.P.RTop.x)/N;
                BR_y(end+1) = OP.P.RTop.y + Y*(OP.P.RBot.y - OP.P.RTop.y)/N;
                if OP.Tols.interior == 1
                    OP.P.Rs.x(end+1) = fR_x(BR_y(end),x_max,PPP,m) - (1/(OP.Tols.sideE+1))*OP.xoff;
                    OP.P.Rs.y(end+1) = BR_y(end);
                end
            end
            BR_ID(1:length(BR_x)) = max(BR_ID) + [1:length(BR_x)];
            BR_x(end+1) = OP.P.RBot.x; BR_y(end+1) = OP.P.RBot.y; BR_ID(end+1) = max(BR_ID) + 1;
            BR_x(end+1) = 0.5*(max(OP.P.Bots.x) + OP.P.RBot.x); BR_y(end+1) = OP.P.RBot.y; BR_ID(end+1) = max(BR_ID) + 1; bts(1:3) = [BR_ID(end-1:end) OP.P.Bots.ID(end)]; bts = fliplr(bts);
            BR_x(end+1:end+length(OP.P.Rights.ID)) = OP.P.Rights.x; BR_y(end+1:end+length(OP.P.Rights.ID)) = OP.P.Rights.y; BR_ID(end+1:end+length(OP.P.Rights.ID)) = OP.P.Rights.ID;
            BR_x(end+1) = 0.5*(max(OP.P.Tops.x) + OP.P.RTop.x); BR_y(end+1) = OP.P.RTop.y; BR_ID(end+1) = max(BR_ID) + 1; tps(1:3) = [BR_ID(1) BR_ID(end) OP.P.Tops.ID(end)]; tps = fliplr(tps);
            [BR_ID, ids] = EXTERNAL_unique_no_sort(BR_ID); BR_x = BR_x(ids); BR_y = BR_y(ids);
            if OP.Tols.interior == 1
                IR_x = OP.P.Rs.x; IR_y = OP.P.Rs.y; IR_ID = [];
                IR_ID(1:length(IR_x)) = max(BR_ID) + [1:length(IR_x)];
            else IR_x = []; IR_y = []; IR_ID = []; end
            
            [TRI_R,VER_R,TRI_R_boundary,VER_R_boundary] = AJSABQ_Delaunay(BR_x,BR_y,BR_ID,IR_x,IR_y,IR_ID,OP.Tols.tolRND,OP.Tols.tolNN,OP.Tols.tolBG,GEO.helix.order,Nmax*2+1,Emax*2+1,OP.Tols.scale);
            disp('Consolidating');
            if GEO.helix.order == 2
                ns = []; bt = find(VER_R.y < min(VER_R.y) + OP.Tols.tolNN); tp = find(VER_R.y > max(VER_R.y) - OP.Tols.tolNN); 
                bts(end+1:end+length(bt)) = VER_R.ID(bt); tps(end+1:end+length(tp)) = VER_R.ID(tp);
                bts = unique(bts); tps = unique(tps); bts_x = []; tps_x = [];
                for bt = 1:length(bts)
                    bts_x(bt) = VER_R.x(find(VER_R.ID == bts(bt)));
                    tps_x(bt) = VER_R.x(find(VER_R.ID == tps(bt)));
                end
                [bts_x, b] = sort(bts_x); [tps_x, t] = sort(tps_x);
                bts = bts(b); tps = tps(t);
            end
            for T = 1:length(TRI_R.ID);
                n1 = find(VER_R.ID == TRI_R.vertices(T,1)); x1 = VER_R.x(n1); y1 = VER_R.y(n1);
                n2 = find(VER_R.ID == TRI_R.vertices(T,2)); x2 = VER_R.x(n2); y2 = VER_R.y(n2);
                n3 = find(VER_R.ID == TRI_R.vertices(T,3)); x3 = VER_R.x(n3); y3 = VER_R.y(n3);
                if GEO.helix.order == 2; WH = [];
                    n4 = find(VER_R.ID == TRI_R.vertices(T,4)); x4 = VER_R.x(n4); y4 = VER_R.y(n4); b4 = 1; if isempty(find(VER_R_boundary == VER_R.ID(n4), 1)); b4 = 0; end
                    n5 = find(VER_R.ID == TRI_R.vertices(T,5)); x5 = VER_R.x(n5); y5 = VER_R.y(n5); b5 = 1; if isempty(find(VER_R_boundary == VER_R.ID(n5), 1)); b5 = 0; end
                    n6 = find(VER_R.ID == TRI_R.vertices(T,6)); x6 = VER_R.x(n6); y6 = VER_R.y(n6); b6 = 1; if isempty(find(VER_R_boundary == VER_R.ID(n6), 1)); b6 = 0; end
                    if b4 == 1; WH(end+1) = 12; end
                    if b5 == 1; WH(end+1) = 23; end
                    if b6 == 1; WH(end+1) = 13; end
                    for W = 1:length(WH)
                        if WH(W) == 12; x_e = x4; y_e = y4; n = n4; N = 4; end
                        if WH(W) == 13; x_e = x6; y_e = y6; n = n6; N = 6; end
                        if WH(W) == 23; x_e = x5; y_e = y5; n = n5; N = 5; end
                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Rights2.x; OP.P.Rights2.y]);
                        PTx = OP.P.Rights2.x(IDX); PTy = OP.P.Rights2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_R.vertices(T,N) = OP.P.Rights2.ID(IDX); ns(end+1) = n; end
                    end
                end
                x_edge12 = 0.5*(x1 + x2); y_edge12 = 0.5*(y1 + y2);
                x_edge13 = 0.5*(x1 + x3); y_edge13 = 0.5*(y1 + y3);
                x_edge23 = 0.5*(x2 + x3); y_edge23 = 0.5*(y2 + y3);
                if abs(x_edge12 - fR_x(y_edge12,x_max,PPP,m)) < OP.Tols.tolNN
                    NSET.Right(end+1) = VER_R.ID(n1); NSET.Right(end+1) = VER_R.ID(n2); ELSET.Right(end+1) = TRI_R.ID(T);
                    if GEO.helix.order == 1; TRI_R.vertices(T,:) = [VER_R.ID(n2) VER_R.ID(n3) VER_R.ID(n1)]; end
                    if GEO.helix.order == 2; TRI_R.vertices(T,:) = [VER_R.ID(n2) VER_R.ID(n3) VER_R.ID(n1) VER_R.ID(n5) VER_R.ID(n6) VER_R.ID(n4)]; NSET.Right(end+1) = VER_R.ID(n4); end
                end
                if abs(x_edge13 - fR_x(y_edge13,x_max,PPP,m)) < OP.Tols.tolNN
                    NSET.Right(end+1) = VER_R.ID(n1); NSET.Right(end+1) = VER_R.ID(n3); ELSET.Right(end+1) = TRI_R.ID(T);
                    if GEO.helix.order == 2; NSET.Right(end+1) = VER_R.ID(n6); end
                end
                if abs(x_edge23 - fR_x(y_edge23,x_max,PPP,m)) < OP.Tols.tolNN
                    NSET.Right(end+1) = VER_R.ID(n2); NSET.Right(end+1) = VER_R.ID(n3); ELSET.Right(end+1) = TRI_R.ID(T);
                    if GEO.helix.order == 1; TRI_R.vertices(T,:) = [VER_R.ID(n3) VER_R.ID(n1) VER_R.ID(n2)]; end
                    if GEO.helix.order == 2; TRI_R.vertices(T,:) = [VER_R.ID(n3) VER_R.ID(n1) VER_R.ID(n2) VER_R.ID(n6) VER_R.ID(n4) VER_R.ID(n5)]; NSET.Right(end+1) = VER_R.ID(n5); end
                end
                for R = 1:length(bts)
                    n = find(TRI_R.vertices(T,:) == bts(R));
                    if ~ isempty(n); TRI_R.vertices(T,n) = tps(R); end
                end
                if (x2 - x1)*(y3 - y1) - (y2 - y1)*(x3 - x1) < 0.0 % Verifying the consistency of the element normal
                   if GEO.helix.order == 1; TRI_R.vertices(T,:) = [TRI_R.vertices(T,3) TRI_R.vertices(T,2) TRI_R.vertices(T,1)]; end
                   if GEO.helix.order == 2; TRI_R.vertices(T,:) = [TRI_R.vertices(T,3) TRI_R.vertices(T,2) TRI_R.vertices(T,1) TRI_R.vertices(T,5) TRI_R.vertices(T,4) TRI_R.vertices(T,6)]; end                   
                end
            end
            if GEO.helix.order == 1; en = 1; end
            if GEO.helix.order == 2
                VER_R.ID(ns) = []; VER_R.x(ns) = []; VER_R.y(ns) = [];
                if GEO.helix.element == 'j'; en = 3; else en = 2; end
            end
            NSET.Top = unique(NSET.Right); NSET = rmfield(NSET,'Right');
            ELSET.Top = ELSET.Right; ELSET = rmfield(ELSET,'Right');

            ELEMENTS.nodes(end+1:end+size(TRI_R.vertices,1),:) = [TRI_R.vertices zeros(size(TRI_R.vertices,1),en)];
            ELEMENTS.index(end+1:end+length(TRI_R.ID)) = TRI_R.ID;
            if GEO.helix.order == 1; ELEMENTS.type(end+1:end+length(TRI_R.ID)) = char(ones(1,length(TRI_R.ID))*'d'); len = length(NODES.index); end
            if GEO.helix.order == 2; ELEMENTS.type(end+1:end+length(TRI_R.ID)) = char(ones(1,length(TRI_R.ID))*'f'); len = length(NODES.index); end
            NODES.index(end+1:end+length(VER_R.ID)) = VER_R.ID;
            NODES.x(end+1:end+length(VER_R.ID)) = VER_R.x;
            NODES.y(end+1:end+length(VER_R.ID)) = VER_R.y;
            for N = 1:length(VER_R.ID)
                % Radial coordinate
                NR = N_R(OP.R,d0,abs(NODES.y(len+N)-0.5*PPP),lambda);
                [NODES.X(len+N), NODES.Y(len+N), NODES.Z(len+N)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODES.x(len+N),NODES.y(len+N),NR);
            end

            MPI_Send(master,2000000,PP.comm,ELEMENTS,NODES,NSET.Top,ELSET.Top);

        elseif isempty(find(slave == rank,1)) == 0 && rank == 2
            [Nmax, Emax, OP, GEO] = MPI_Recv(master,1500000,PP.comm); % Receiving info from master thread about RHS bit
            NODES = struct('index',[],'x',[],'y',[],'X',[],'Y',[],'Z',[],'case',[]);
            ELEMENTS = struct('index',[],'nodes',[],'type',[],'case',[]);

            % CPU >= 3 - Slave 2 Thread - Triangulating the misalingment
            BB_ID = [OP.P.Bots.ID OP.P.Ints.ID]; BB_x = [OP.P.Bots.x OP.P.Ints.x]; BB_y = [OP.P.Bots.y OP.P.Ints.y]; ns = [];

            [TRI_B,VER_B,TRI_B_boundary,VER_B_boundary] = AJSABQ_Delaunay(BB_x,BB_y,BB_ID,[],[],[],OP.Tols.tolRND,OP.Tols.tolNN,OP.Tols.tolBG,GEO.helix.order,Nmax*3+1,Emax*3+1,0);
            disp('Consolidating');
            for T = 1:length(TRI_B.ID)
                if GEO.helix.order == 2; WH = [];
                    n4 = find(VER_B.ID == TRI_B.vertices(T,4)); x4 = VER_B.x(n4); y4 = VER_B.y(n4); b4 = 1; if isempty(find(VER_B_boundary == VER_B.ID(n4), 1)); b4 = 0; end
                    n5 = find(VER_B.ID == TRI_B.vertices(T,5)); x5 = VER_B.x(n5); y5 = VER_B.y(n5); b5 = 1; if isempty(find(VER_B_boundary == VER_B.ID(n5), 1)); b5 = 0; end
                    n6 = find(VER_B.ID == TRI_B.vertices(T,6)); x6 = VER_B.x(n6); y6 = VER_B.y(n6); b6 = 1; if isempty(find(VER_B_boundary == VER_B.ID(n6), 1)); b6 = 0; end
                    if b4 == 1; WH(end+1) = 12; end
                    if b5 == 1; WH(end+1) = 23; end
                    if b6 == 1; WH(end+1) = 13; end
                    for W = 1:length(WH)
                        if WH(W) == 12; x_e = x4; y_e = y4; n = n4; N = 4; x_o1 = x5; y_o1 = y5; n_o1 = n5; N_o1 = 5; x_o2 = x6; y_o2 = y6; n_o2 = n6; N_o2 = 6; end
                        if WH(W) == 13; x_e = x6; y_e = y6; n = n6; N = 6; x_o1 = x4; y_o1 = y4; n_o1 = n4; N_o1 = 4; x_o2 = x5; y_o2 = y5; n_o2 = n5; N_o2 = 5; end
                        if WH(W) == 23; x_e = x5; y_e = y5; n = n5; N = 5; x_o1 = x4; y_o1 = y4; n_o1 = n4; N_o1 = 4; x_o2 = x6; y_o2 = y6; n_o2 = n6; N_o2 = 6; end
                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Ints2.x; OP.P.Ints2.y]);
                        PTx = OP.P.Ints2.x(IDX); PTy = OP.P.Ints2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Ints2.ID(IDX); ns(end+1) = n; end

                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Bots2.x; OP.P.Bots2.y]);
                        PTx = OP.P.Bots2.x(IDX); PTy = OP.P.Bots2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Bots2.ID(IDX); ns(end+1) = n; end

                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Lefts2.x; OP.P.Lefts2.y]);
                        PTx = OP.P.Lefts2.x(IDX); PTy = OP.P.Lefts2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Lefts2.ID(IDX); ns(end+1) = n; end

                        IDX = EXTERNAL_nearest_neighbour([x_e; y_e],[OP.P.Rights2.x; OP.P.Rights2.y]);
                        PTx = OP.P.Rights2.x(IDX); PTy = OP.P.Rights2.y(IDX);
                        dist = sqrt( (PTx - x_e)*(PTx - x_e) + (PTy - y_e)*(PTy - y_e) );
                        if dist < OP.Tols.tolNN; TRI_B.vertices(T,N) = OP.P.Rights2.ID(IDX); ns(end+1) = n; end
                    end
                end
            end
            for T = 1:length(TRI_B.ID)
                for N = 1:length(OP.P.Bots.ID)
                    a = find(TRI_B.vertices(T,:) == OP.P.Bots.ID(N));
                    if ~ isempty(a); TRI_B.vertices(T,a) = OP.P.Tops.ID(N); end
                end
                for N2 = 1:length(OP.P.Bots2.ID)
                    a2 = find(TRI_B.vertices(T,:) == OP.P.Bots2.ID(N2));
                    if ~ isempty(a2); TRI_B.vertices(T,a2) = OP.P.Tops2.ID(N2); end
                end
            end
            if GEO.helix.order == 1; en = 1; end
            if GEO.helix.order == 2; VER_B.ID(ns) = []; VER_B.x(ns) = []; VER_B.y(ns) = [];
                if GEO.helix.element == 'j'; en = 3; else en = 2; end
            end

            ELEMENTS.nodes(end+1:end+size(TRI_B.vertices,1),:) = [TRI_B.vertices zeros(size(TRI_B.vertices,1),en)];
            ELEMENTS.index(end+1:end+length(TRI_B.ID)) = TRI_B.ID;
            if GEO.helix.order == 1; ELEMENTS.type(end+1:end+length(TRI_B.ID)) = char(ones(1,length(TRI_B.ID))*'d'); len = length(NODES.index); end
            if GEO.helix.order == 2; ELEMENTS.type(end+1:end+length(TRI_B.ID)) = char(ones(1,length(TRI_B.ID))*'f'); len = length(NODES.index); end
            NODES.index(end+1:end+length(VER_B.ID)) = VER_B.ID;
            NODES.x(end+1:end+length(VER_B.ID)) = VER_B.x;
            NODES.y(end+1:end+length(VER_B.ID)) = VER_B.y;
            for N = 1:length(VER_B.ID)
                % Radial coordinate
                NR = N_R(OP.R,d0,abs(NODES.y(len+N)-0.5*PPP),lambda);
                [NODES.X(len+N), NODES.Y(len+N), NODES.Z(len+N)] = AJSABQ_helical_mapping(GEO.helix,GEO.feature,OP.Hc,OP.R,NODES.x(len+N),NODES.y(len+N),NR);
            end
            
            MPI_Send(master,2500000,PP.comm,ELEMENTS,NODES); NSET = []; ELSET = [];
        end

        % Calibrating
        if isempty(find(master == rank,1)) == 0
            % From slave 1 thread - RHS region
            [SlaveELEMENTS, SlaveNODES, SlaveNTop, SlaveETop] = MPI_Recv(slave(1),2000000,PP.comm);
            ELEMENTS.nodes(end+1:end+length(SlaveELEMENTS.index),:) = SlaveELEMENTS.nodes;
            ELEMENTS.index(end+1:end+length(SlaveELEMENTS.index)) = SlaveELEMENTS.index;
            NODES.index(end+1:end+length(SlaveNODES.index)) = SlaveNODES.index;
            NODES.x(end+1:end+length(SlaveNODES.x)) = SlaveNODES.x;
            NODES.y(end+1:end+length(SlaveNODES.y)) = SlaveNODES.y;
            NODES.X(end+1:end+length(SlaveNODES.X)) = SlaveNODES.X;
            NODES.Y(end+1:end+length(SlaveNODES.Y)) = SlaveNODES.Y;
            NODES.Z(end+1:end+length(SlaveNODES.Z)) = SlaveNODES.Z;
            NSET.Top = SlaveNTop; ELSET.Top = SlaveETop;

            % From slave 2 thread - misalignment
            [SlaveELEMENTS, SlaveNODES] = MPI_Recv(slave(2),2500000,PP.comm);
            ELEMENTS.nodes(end+1:end+length(SlaveELEMENTS.index),:) = SlaveELEMENTS.nodes;
            ELEMENTS.index(end+1:end+length(SlaveELEMENTS.index)) = SlaveELEMENTS.index;
            NODES.index(end+1:end+length(SlaveNODES.index)) = SlaveNODES.index;
            NODES.x(end+1:end+length(SlaveNODES.x)) = SlaveNODES.x;
            NODES.y(end+1:end+length(SlaveNODES.y)) = SlaveNODES.y;
            NODES.X(end+1:end+length(SlaveNODES.X)) = SlaveNODES.X;
            NODES.Y(end+1:end+length(SlaveNODES.Y)) = SlaveNODES.Y;
            NODES.Z(end+1:end+length(SlaveNODES.Z)) = SlaveNODES.Z;
            OP.Hel_Tri = length(ELEMENTS.index) - OP.Hel_Quad;
        end
end