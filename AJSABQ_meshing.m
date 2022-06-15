function [NODES,ELEMENTS,NSET,OP] = AJSABQ_meshing(PP,RSL,OP,NODES,ELEMENTS,NSET,GEO,MATERIAL)
% MPI Parallelised meshing algorithm
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:24 (previously 07/08/11 - 19:40)

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

siz = []; RemT = []; RemTemp = []; RemTOT = []; p1 = pi/180; disp(['Meshing (shell elements): Rank ',num2str(PP.rank)]);
i1 = OP.i1; i2 = OP.i2; e1 = OP.e1; e2 = OP.e2; div = OP.div; Hmin = OP.Hmin; Hmax = OP.Hmax; Ttot = OP.T;
master = PP.master; slave = PP.slave; rank = PP.rank; comm = PP.comm; CPUs = PP.CPUs; tp = 0;

for G = 1:length(RSL.Region)
    for D = 1:CPUs
        if isempty(find(master == rank,1)) == 0; D = 1; end
        if isempty(find(slave == rank,1)) == 0; D = PP.rank + 1;
            NODES = struct('index',[],'R',[],'T',[],'Z',[]);
            ELEMENTS = struct('index',[],'nodes',[],'region',[],'type',[],'t',[],'alt',[]);
        end
    end
    siz(length(siz)+1) = length(RemT); RemT = []; RTcurrent = RSL.Ttot+1; alpres = 0;
    Imin = i1(G,D); Imax = i2(G,D); hmin = Hmin(G,D); Emin = e1(G,D); Emax = e2(G,D); hmax = Hmax(G,D); I = Imin - 1; E = Emin - 1;
    if G == 1 && D == 1; Rno = abs(hmax - hmin)/div(G,D); else Rno = abs(hmax - hmin)/(div(G,D)-1); end

    for Z = hmin:Rno:hmax
        for A = 1:length(RSL.t1)
            tmax = RSL.t2(A); tmin = RSL.t1(A); Tno = abs(tmax - tmin)/RSL.T(A); ZLalt = 0;
            if A == 1; tmin = tmin; else tmin = tmin + Tno; end
            for T = tmin:Tno:tmax; I = I + 1;
                if (Ttot == 360 && T == 360 && isempty(find(master == rank,1)) == 0); RemT(length(RemT)+1) = I;
                    NODES.index(I) = I; NODES.R(I) = 0; NODES.T(I) = 0; NODES.Z(I) = 0; continue
                end
                if (Ttot == 360 && T == 360 && isempty(find(slave == rank,1)) == 0); RemT(length(RemT)+1) = I;
                    NODES.index(I+1-Imin) = I; NODES.R(I+1-Imin) = 0; NODES.T(I+1-Imin) = 0; NODES.Z(I+1-Imin) = 0; continue
                end

                %%% PRIMARY GEOMETRIC FEATURE (i.e. cone/cylinder)
                % Hopper (if present)
                if GEO.hop.Inc == 1 && RSL.Where(G) == 1; NR = (Z/OP.Hh)*OP.R; end

                % Cylinder (if present)
                if GEO.cyl.Inc == 1 && RSL.Where(G) == 2; NR = OP.R; end

                % Roof (if present)
                if GEO.roof.Inc == 1 && RSL.Where(G) == 3; NR = OP.R - (Z - OP.H)/tan(GEO.roof.Angle*pi/180); end

                %%% SECONDARY GEOMETRIC FEATURE (i.e. weld/perturbation/indentation/flattening/hopper plate)
                % Circumferential Welds (if present)
                if GEO.feature.C_Weld.Toggle == 1
                    delta = pi*sqrt(OP.R*RSL.t(G))/(3*(1-MATERIAL.v*MATERIAL.v))^0.25;
                    for W = 1:length(GEO.feature.C_Weld.Type); ZR = abs(GEO.feature.C_Weld.Z(W) - Z);
                        if GEO.feature.C_Weld.Type(W) == 'A'; s = 1; elseif GEO.feature.C_Weld.Type(W) == 'B'; s = 0; end
                        if G == length(RSL.Region); ampt = GEO.feature.C_Weld.Ampt(W)*RSL.t(G); else ampt = GEO.feature.C_Weld.Ampt(W)*RSL.t(G+1); end
                        kkk = 10^(-1); nnn = 10;
                        NR = NR - ampt*exp(-ZR*pi/delta)*(cos(pi*ZR/delta)+s*sin(pi*ZR/delta))*(1+kkk*cos(nnn*T*pi/180));
                    end
                end

                % Axial Welds (if present)
                if GEO.feature.A_Weld.Toggle == 1
                    for ZL = 1:length(GEO.feature.A_Weld.Z) - 1
                        if mod(ZL,2) == 1; ZLalt = 1; elseif mod(ZL,2) == 0; ZLalt = 2; end
                        if GEO.feature.A_Weld.Type == 'A'; s = 1; elseif GEO.feature.A_Weld.Type == 'B'; s = 0; end
                        deltaT = (pi*sqrt(OP.R*RSL.t(G))/(3*(1-MATERIAL.v*MATERIAL.v))^0.25)*pi/180;
                        ampT = GEO.feature.A_Weld.Ampt(find(GEO.feature.A_Weld.t == RSL.t(G)))*RSL.t(G);
                        if Z >= GEO.feature.A_Weld.Z(ZL) && Z <= GEO.feature.A_Weld.Z(ZL+1)
                            for TT = 0:GEO.feature.A_Weld.Tfreq:OP.T
                                if ZLalt == 1; TR = abs(TT - T); elseif ZLalt == 2; TR = abs(TT + GEO.feature.A_Weld.Tphase - T); end
                                NR = NR - ampT*exp(-TR*pi/deltaT)*(cos(pi*TR/deltaT)+s*sin(pi*TR/deltaT));
                            end
                        end
                    end
                end

                % Super-elliptical flattening (if present)
                if GEO.feature.Super_ell(1) == 1 && RSL.Where(G) == 2
                    e_d = GEO.feature.Super_ell(2); tS = GEO.feature.Super_ell(3);
                    e_p = GEO.feature.Super_ell(4); e_q = GEO.feature.Super_ell(5); Zo = GEO.feature.Super_ell(6);
                    if T <= tS || T >= (360 - tS); Tell = (T/tS)*90;
                        if GEO.feature.Super_ell(7) == 1; kk = 1; else kk = 0; end
                        e_a = OP.R*sin(tS*pi/180); e_b = OP.R*(1 - cos(tS*pi/180)) - e_d; save('e_a','e_b');
                        e_x = e_a*sin(Tell*pi/180)^(2/e_p); e_y = e_b*cos(Tell*pi/180)^(2/e_q) + OP.R*cos(tS*pi/180);
                        dep = OP.R - sqrt(e_x*e_x + e_y*e_y);
                        if Z - OP.Hh <= Zo; NR = NR - dep*sin(((Z - OP.Hh)/(2*Zo))*pi)*kk;
                        elseif Z - OP.Hh >= Zo; NR = NR - dep*sin(((OP.Hc - Z - OP.Hh)/(2*(OP.Hc - Zo)))*pi)*kk;
                        end
                    end
                end             
                
                % Circular flattening (if present)
                if GEO.feature.Circ_flat(1) == 1 && RSL.Where(G) == 2
                    d = GEO.feature.Circ_flat(2); tC = GEO.feature.Circ_flat(3)*pi/180; dep = 0;
                    Zo = GEO.feature.Circ_flat(4);
                    DD = -(1/2)*d*(2*OP.R-d)/(-OP.R+d+OP.R*cos(tC));
                    Rp = (1/2)*(-2*OP.R^2+2*OP.R*d-d^2+2*cos(tC)*OP.R^2-2*cos(tC)*OP.R*d)/(-OP.R+d+OP.R*cos(tC));
                    tCP = (180/pi)*atan2(2*OP.R*sin(tC)*(-OP.R+d+OP.R*cos(tC))/(-2*OP.R^2+2*OP.R*d-d^2+2*cos(tC)*OP.R^2-2*cos(tC)*OP.R*d), (-2*OP.R*d+2*OP.R^2*cos(tC)^2+d^2+2*cos(tC)*OP.R*d-2*cos(tC)*OP.R^2)/(-2*OP.R^2+2*OP.R*d-d^2+2*cos(tC)*OP.R^2-2*cos(tC)*OP.R*d));
                    if T <= tC*180/pi || T >= (360 - tC*180/pi)
                        Tfrac = T/(tC*180/pi); tP = Tfrac*tCP;
                        x = Rp*cos(tP*pi/180) - DD; y = Rp*sin(tP*pi/180); dep = OP.R - sqrt(x*x+y*y);
                    end
                    if GEO.feature.Circ_flat(5) == 1; kk = 1; else kk = 0; end
                    if Z - OP.Hh <= Zo; NR = NR - dep*sin(((Z - OP.Hh)/(2*Zo))*pi)*kk;
                    elseif Z - OP.Hh >= Zo; NR = NR - dep*sin(((OP.Hc - Z - OP.Hh)/(2*(OP.Hc - Zo)))*pi)*kk;
                    end
                end

                % Inward sine wave (if present)
                if GEO.feature.Axi_sin(1) == 1 && RSL.Where(G) == 2; NR = NR - GEO.feature.Axi_sin(2)*sin(((Z - OP.Hh)/OP.Hc)*pi*GEO.feature.Axi_sin(3)); end

                % Inward circumferential cosine wave (if present)
                if GEO.feature.Circ_cos(1) == 1 && RSL.Where(G) == 2
                    %if GEO.feature.Circ_cos(4) == 1; kk = sin((0.5*(Z - OP.Hh)/OP.Hc)*pi); else kk = 1; end
                    if GEO.feature.Circ_cos(4) == 1; kk = sin((0.5*(OP.Hc - Z)/OP.Hc)*pi); else kk = 1; end
                    NR = NR - GEO.feature.Circ_cos(2)*cos(GEO.feature.Circ_cos(3)*pi*T/180)*kk;
                end

                % Eigenmode mesh perturbation (if present)
                % HACKY - IT IS 2022 NOW AND I DON'T REMEMBER WHAT THIS DOES
                if GEO.feature.Cyl_perturbation(1) == 1 && RSL.Where(G) == 2
                    %NR = NR - GEO.feature.Cyl_perturbation(4)*sin(GEO.feature.Cyl_perturbation(2)*T*pi/180)*sin(0.5*GEO.feature.Cyl_perturbation(3)*((OP.Hc - Z)/OP.Hc)*(2*pi));
                    %
                    %                     % Temporary alternative
                    %                     Fourier_Zs = [-0.001059438, -0.745805852, 0.002496493, -0.508437432, 1.133635735, -0.31955094, -0.924383518, -0.222963867, -0.126870896, -0.126689242, 0.016646482, -0.054526238,...
                    %                         0.049798784, -0.010960989, 0.043399938, 0.008505051, 0.025813068, 0.012291749, 0.010344266, 0.008917276, 0.001405561, 0.004364052, -0.001739362, 0.001343646, -0.001723206,...
                    %                         0.000119111, -0.000849606, -9.34857E-05, -0.000236359, -3.23756E-05, -2.59669E-05, 1.54725E-06, 4.76677E-07, -3.15937E-07, 1.56806E-07, 8.9212E-08, -3.91361E-07, 1.15409E-07,...
                    %                         1.4054E-06, -7.48406E-06, -1.39563E-05, -9.75194E-05, -6.02569E-05, -0.000411019, -1.65092E-05, -0.000961072, 0.000502677, -0.001267055, 0.00201881, -0.000120809, 0.004646746,...
                    %                         0.00412271, 0.007268782, 0.012460928, 0.006865736, 0.023449624, -0.001404409, 0.031005494, -0.023080593, 0.021181368, -0.06260611, -0.039647505, -0.121797201, -0.334795892,...
                    %                         -0.207696364, 1.006396786, -0.191452544, 0.137497333, -0.36360562, 0.04707476, 0.026574979, 0.312203168, -0.467739967, -0.437036117, -0.339123156, 1.13036779, -0.081616182,...
                    %                         0.110252563, -0.125352113, 0.005243791, -0.079492413, -0.023677334, -0.038730378, -0.024858204, -0.011562658, -0.016530341, 0.002205812, -0.007663445, 0.006290986, -0.001869571,...
                    %                         0.005304738, 0.000583571, 0.002910503, 0.00093786, 0.00105195, 0.000542389, 0.000179099, 0.000175514, -3.53462E-05, 2.52552E-05, -2.20767E-05]';
                    %
                    %                     rad = 0; L = length(Fourier_Zs); INTS = [0]; c = 0;
                    %                     for l = 1:L/2; c = c + 1;
                    %                         INTS(end+1:end+2) = c*[1 1];
                    %                     end
                    %                     if length(INTS) > L; INTS(end) = []; end
                    %                     if length(INTS) >= 1; coss = [1 2:2:length(INTS)]; else coss = [1]; end
                    %                     if length(INTS) >= 3; sins = [3:2:length(INTS)]; else sins = []; end
                    %
                    %                     rad = rad + sum(Fourier_Zs(coss).*cos(INTS(coss)'*Z*pi/180));
                    %                     rad = rad + sum(Fourier_Zs(sins).*sin(INTS(sins)'*Z*pi/180));
                    load AJS_spline
                    rad = fittedmodel1(Z);
                    if Z <= 1e-3 || Z >= OP.Hc - 1e-3; rad = 0; end
                    NR = NR+GEO.feature.Cyl_perturbation(4)*sin(GEO.feature.Cyl_perturbation(2)*T*pi/180)*rad;
                end

                %%% ASSIGNING NODES TO SETS
                if isempty(find(master == rank,1)) == 0; NODES.index(I) = I; NODES.R(I) = NR; NODES.T(I) = T; NODES.Z(I) = Z; end
                if isempty(find(slave == rank,1)) == 0; NODES.index(I+1-Imin) = I; NODES.R(I+1-Imin) = NR; NODES.T(I+1-Imin) = T; NODES.Z(I+1-Imin) = Z; end
                if Z == OP.Hh; NSET.Bottom(length(NSET.Bottom)+1) = I; end
                if GEO.hop.Inc == 1 && Z == min(RSL.h1); NSET.Czelusc(length(NSET.Czelusc)+1) = I; end
                if GEO.cyl.Inc == 1 && Z == hmax; NSET.Top(length(NSET.Top)+1) = I; end
                if T ~= 360 && T == 0; NSET.Right(length(NSET.Right)+1) = I; end
                if T ~= 360 && T == OP.T; NSET.Left(length(NSET.Left)+1) = I; end
            end
        end
    end

    if Ttot == 360
        correction = RemT - RSL.Ttot; RemTemp = RemT;
        if G > 1 && isempty(find(master == rank,1)) == 0; RemTemp(length(RemT) + 1) = RemTemp(1) - RSL.Ttot - 1; correction(length(correction) + 1) = RemTemp(1) - 2*RSL.Ttot - 1; end
        if isempty(find(slave == rank,1)) == 0; RemTemp(length(RemTemp) + 1) = RemTemp(1) - RSL.Ttot - 1; correction(length(correction) + 1) = RemTemp(1) - 2*RSL.Ttot - 1; end
    end
    if (GEO.roof.Inc == 1 && RSL.Where(G) == 3 && CPUs == 1) || (GEO.roof.Inc == 1 && RSL.Where(G) == 3 && CPUs > 1 && rank == max(slave))
        RemR = NODES.index((length(NODES.index)-RSL.Ttot+tp):length(NODES.index));
        if Ttot == 360
            RemTemp(length(RemTemp)+1:length(RemTemp)+length(RemR)) = RemR;
            correction(length(correction)+1:length(correction)+length(RemR)) = (max(NODES.index)-1)*ones(1,length(RemR));
            RemT(length(RemT)+1:length(RemT)+length(RemR)) = RemR;
        else RemT = RemR; RemTemp = RemR;
            correction = max(NODES.index)*ones(1,length(RemR));
        end
        if Ttot == 360; RemTemp(length(RemTemp)) = []; RemT(length(RemT)) = []; correction(length(correction)) = []; end
    end
    if G == 1 && CPUs == 1; BaseZ = 1; elseif G > 1 && CPUs == 1; BaseZ = sum(RSL.Z(1:G-1)) + 1; end
    if CPUs == 1; TopZ1 = BaseZ + RSL.Z(G) - 1; TopZ2 = BaseZ + RSL.Z(G) - 2; end
    if G == 1 && CPUs > 1; BaseZ = 1 + sum(div(G,1:rank)); elseif G > 1 && CPUs > 1; BaseZ = sum(RSL.Z(1:G-1)) + sum(div(G,1:rank)) + 1; end
    if CPUs > 1; TopZ1 = BaseZ + div(G,rank+1) - 1; TopZ2 = BaseZ + div(G,rank+1) - 2; end
    if isempty(find(master == rank,1)) == 0; J = E; end
    if isempty(find(slave == rank,1)) == 0; J = E + 1 - Emin; end

    %%% ELEMENT GENERATION
    if RSL.Type(G) == 'a' || RSL.Type(G) == 'b' || RSL.Type(G) == 'c' % If it is a 4-node element
        for Zn = BaseZ:TopZ1 % For every row but one of vertical nodes
            for Tn = 1:RSL.Ttot; E = E + 1; J = J + 1; % For every column but one of circumferential nodes
                ELEMENTS.region(J) = RSL.Region(G); ELEMENTS.type(J) = RSL.Type(G); ELEMENTS.index(J) = E;
                ELEMENTS.nodes(J,1) = Tn + (Zn-1)*RTcurrent; ELEMENTS.nodes(J,2) = Tn + 1 + (Zn-1)*RTcurrent;
                ELEMENTS.nodes(J,3) = Tn + 1 + Zn*RTcurrent; ELEMENTS.nodes(J,4) = Tn + Zn*RTcurrent;
                ELEMENTS.t(J) = RSL.t(G);
                if (Ttot == 360 || (GEO.roof.Inc == 1 && RSL.Where(G) == 3))
                    for K = 1:length(RemTemp)
                        clear('A'); A = find(ELEMENTS.nodes(J,:) == RemTemp(K)); if isempty(A) ~= 1; ELEMENTS.nodes(J,A) = correction(K); end
                    end
                end
            end
        end

    elseif RSL.Type(G) == 'd' || RSL.Type(G) == 'e' % If it is a 3-node element
        for Zn = BaseZ:TopZ1 % For every row but one of vertical nodes
            for Tn = 1:RSL.Ttot; E = E + 1; J = J + 1; % For every column but one of circumferential nodes
                % Is it the triangular element with 2 nodes at bottom (alt = 1) or 2 nodes at top (alt = 2)
                ELEMENTS.region(J) = RSL.Region(G); ELEMENTS.type(J) = RSL.Type(G); ELEMENTS.index(J) = E; ELEMENTS.alt(J) = 1; ELEMENTS.t(J) = RSL.t(G);
                ELEMENTS.nodes(J,1) = Tn + (Zn-1)*RTcurrent; ELEMENTS.nodes(J,2) = Tn + 1 + (Zn-1)*RTcurrent; ELEMENTS.nodes(J,3) = Tn + 1 + Zn*RTcurrent;
                if (Ttot == 360 || (GEO.roof.Inc == 1 && RSL.Where(G) == 3))
                    for K = 1:length(RemTemp)
                        clear('A'); A = find(ELEMENTS.nodes(J,:) == RemTemp(K)); if isempty(A) ~= 1; ELEMENTS.nodes(J,A) = correction(K); end
                    end
                end
                E = E + 1; J = J + 1; ELEMENTS.region(J) = RSL.Region(G); ELEMENTS.type(J) = RSL.Type(G); ELEMENTS.index(J) = E; ELEMENTS.alt(J) = 2;
                ELEMENTS.nodes(J,1) = Tn + 1 + Zn*RTcurrent; ELEMENTS.nodes(J,2) = Tn + Zn*RTcurrent; ELEMENTS.nodes(J,3) = Tn + (Zn-1)*RTcurrent; ELEMENTS.t(J) = RSL.t(G);
                if (Ttot == 360 || (GEO.roof.Inc == 1 && RSL.Where(G) == 3))
                    for K = 1:length(RemTemp)
                        clear('A'); A = find(ELEMENTS.nodes(J,:) == RemTemp(K)); if isempty(A) ~= 1; ELEMENTS.nodes(J,A) = correction(K); end
                    end
                end
            end
        end

    elseif RSL.Type(G) == 'f' % If it is a 6-node element
        for Zn = BaseZ:2:TopZ2
            for Tn = 1:2:RSL.Ttot; E = E + 1; J = J + 1;
                ELEMENTS.region(J) = RSL.Region(G); ELEMENTS.type(J) = RSL.Type(G); ELEMENTS.index(J) = E; ELEMENTS.alt(J) = 1; ELEMENTS.t(J) = RSL.t(G);
                ELEMENTS.nodes(J,1) = Tn + (Zn-1)*RTcurrent; ELEMENTS.nodes(J,2) = Tn + 2 + (Zn-1)*RTcurrent;
                ELEMENTS.nodes(J,3) = Tn + 2 + (Zn+1)*RTcurrent; ELEMENTS.nodes(J,4) = Tn + 1 + (Zn-1)*RTcurrent;
                ELEMENTS.nodes(J,5) = Tn + 2 + Zn*RTcurrent; ELEMENTS.nodes(J,6) = Tn + 1 + Zn*RTcurrent;
                if (Ttot == 360 || (GEO.roof.Inc == 1 && RSL.Where(G) == 3))
                    for K = 1:length(RemTemp)
                        clear('A'); A = find(ELEMENTS.nodes(J,:) == RemTemp(K)); if isempty(A) ~= 1; ELEMENTS.nodes(J,A) = correction(K); end
                    end
                end
                E = E + 1; J = J + 1; ELEMENTS.region(J) = RSL.Region(G); ELEMENTS.type(J) = RSL.Type(G); ELEMENTS.index(J) = E; ELEMENTS.alt(J) = 2;
                ELEMENTS.nodes(J,1) = Tn + 2 + (Zn+1)*RTcurrent; ELEMENTS.nodes(J,2) = Tn + (Zn+1)*RTcurrent; ELEMENTS.t(J) = RSL.t(G);
                ELEMENTS.nodes(J,3) = Tn + (Zn-1)*RTcurrent; ELEMENTS.nodes(J,4) = Tn + 1 + (Zn+1)*RTcurrent;
                ELEMENTS.nodes(J,5) = Tn + Zn*RTcurrent; ELEMENTS.nodes(J,6) = Tn + 1 + Zn*RTcurrent;
                if (Ttot == 360 || (GEO.roof.Inc == 1 && RSL.Where(G) == 3))
                    for K = 1:length(RemTemp)
                        clear('A'); A = find(ELEMENTS.nodes(J,:) == RemTemp(K)); if isempty(A) ~= 1; ELEMENTS.nodes(J,A) = correction(K); end
                    end
                end
            end
        end

    elseif RSL.Type(G) == 'g' || RSL.Type(G) == 'h' || RSL.Type(G) == 'j' % If it is an 8-node or 9-node element
        for Zn = BaseZ:2:TopZ2 %(BaseZ + RSL.Z(G)-2); % For every row but one of vertical nodes
            for Tn = 1:2:RSL.Ttot; E = E + 1; J = J + 1; % For every column but one of circumferential nodes
                ELEMENTS.region(J) = RSL.Region(G); ELEMENTS.type(J) = RSL.Type(G); ELEMENTS.index(J) = E;
                ELEMENTS.nodes(J,1) = Tn + (Zn-1)*RTcurrent; ELEMENTS.nodes(J,2) = Tn + 2 + (Zn-1)*RTcurrent;
                ELEMENTS.nodes(J,3) = Tn + 2 + (Zn+1)*RTcurrent; ELEMENTS.nodes(J,4) = Tn + (Zn+1)*RTcurrent;
                ELEMENTS.nodes(J,5) = Tn + 1 + (Zn-1)*RTcurrent; ELEMENTS.nodes(J,6) = Tn + 2 + Zn*RTcurrent;
                ELEMENTS.nodes(J,7) = Tn + 1 + (Zn+1)*RTcurrent; ELEMENTS.nodes(J,8) = Tn + Zn*RTcurrent;
                if RSL.Type(G) == 'j'; ELEMENTS.nodes(J,9) = Tn + 1 + Zn*RTcurrent; end
                ELEMENTS.t(J) = RSL.t(G);
                if (Ttot == 360 || (GEO.roof.Inc == 1 && RSL.Where(G) == 3))
                    for K = 1:length(RemTemp)
                        clear('A'); A = find(ELEMENTS.nodes(J,:) == RemTemp(K)); if isempty(A) ~= 1; ELEMENTS.nodes(J,A) = correction(K); end
                    end
                end
            end
        end
    end

    %%% CALIBRATION
    if isempty(find(slave == rank,1)) == 0; MPI_Send(master,4+G+rank,comm,NODES,ELEMENTS,RemTemp,NSET); NSET = struct('Bottom',[],'Top',[],'Right',[],'Left',[],'Czelusc',[]); end
    if isempty(find(master == rank,1)) == 0 && CPUs > 1
        for S = 1:length(slave)
            [SlaveNodes, SlaveElements, RemDT, SlaveNSet] = MPI_Recv(slave(S),4+G+slave(S),comm);
            ELEMENTS.index(SlaveElements.index) = SlaveElements.index; ELEMENTS.t(SlaveElements.index) = SlaveElements.t;
            for L = 1:size(SlaveElements.nodes,2)
                ELEMENTS.nodes(SlaveElements.index,L) = SlaveElements.nodes(:,L);
            end
            ELEMENTS.region(SlaveElements.index) = SlaveElements.region; ELEMENTS.type(SlaveElements.index) = SlaveElements.type;
            if RSL.Type(G) == 'd' || RSL.Type(G) == 'e' || RSL.Type(G) == 'f'; ELEMENTS.alt(SlaveElements.index) = SlaveElements.alt; end
            NODES.index(SlaveNodes.index) = SlaveNodes.index; NODES.Z(SlaveNodes.index) = SlaveNodes.Z;
            NODES.T(SlaveNodes.index) = SlaveNodes.T; NODES.R(SlaveNodes.index) = SlaveNodes.R;
            if ~isempty(SlaveNSet.Bottom > 0); NSET.Bottom(length(NSET.Bottom)+1:length(NSET.Bottom)+length(SlaveNSet.Bottom)) = SlaveNSet.Bottom; end
            if ~isempty(SlaveNSet.Czelusc > 0); NSET.Czelusc(length(NSET.Czelusc)+1:length(NSET.Czelusc)+length(SlaveNSet.Czelusc)) = SlaveNSet.Czelusc; end
            if ~isempty(SlaveNSet.Top > 0); NSET.Top(length(NSET.Top)+1:length(NSET.Top)+length(SlaveNSet.Top)) = SlaveNSet.Top; end
            if Ttot ~= 360 && isempty(NSET.Right) == 0; NSET.Right(length(NSET.Right)+1:length(NSET.Right)+length(SlaveNSet.Right)) = SlaveNSet.Right; end
            if Ttot ~= 360 && isempty(NSET.Left) == 0; NSET.Left(length(NSET.Left)+1:length(NSET.Left)+length(SlaveNSet.Left)) = SlaveNSet.Left; end
            RemT(length(RemT)+1:length(RemT)+length(RemDT)) = RemDT;
            RemTOT(length(RemTOT)+1:length(RemTOT)+length(RemT)) = RemT;
        end
    end
    if isempty(find(master == rank,1)) == 0 && CPUs == 1
        RemTOT(length(RemTOT)+1:length(RemTOT)+length(RemTemp)) = RemTemp;
    end
end
OP.RemTOT = RemTOT;