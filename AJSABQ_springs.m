function [NODES,SPRING,Axis] = AJSABQ_springs(NODES,SOLID,SPRING,OP,Axis,TYP,Ttot)
% EXPERIMENTAL spring element meshing algorithm 
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:27 (previously 28/01/10 - 00:06)

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

%SPRING = struct('toggle',[],'index',[],'nodes',[],'region',[],'set',[],'constant',[],'discharge',[]); % Struct listing spring elements

disp('Meshing (spring elements).');
Ntemp = struct('index',[],'R',[],'T',[],'Z',[]);

Noffset = max(NODES.index);
Eoffset = OP.Hop + OP.Cyl + OP.Roof;
n = 0; nS = 0; J = 0; S = 0; SET = 0;

SPRING.setlist = struct('no',[],'k',[]);

% Calculating values necessary for local pressures
if SPRING.Discharge == 1
    if TYP == 19
        Frang = SOLID.Frang*pi/180; SOLID.Flowr = OP.kc(1)*OP.R; yk = tan(Frang)*SOLID.K/(SOLID.mewL*SOLID.K);
        vals = [OP.kc(1), yk]; tol = 0.01; N = 4; O = 0; P = 0; ex = 0; TH = []; PS = [];
        while ex == 0; O = O + 1;
            if O <= 4; go = 0;
                while go ~= N; go = 0;
                    inp1 = rand*pi/2; inp2 = rand*pi/2;
                    ans = fminsearch(@(guess) AJSABQ_objective_function(guess,vals),[inp1;inp2]);
                    THc = ans(1); Psi = ans(2);
                    ecoR = cos(THc) - OP.kc*cos(Psi);
                    if sin(THc)/sin(Psi) <= OP.kc + tol && sin(THc)/sin(Psi) >= OP.kc - tol; go = go + 1; else go = 0; end
                    if Psi*180/pi < 90 && THc*180/pi < 90; go = go + 1; else go = 0; end
                    if Psi >= THc; go = go + 1; else go = 0; end
                    if Psi >= 0 && THc >= 0; go = go + 1; else go = 0; end
                    if go == N; break; end
                end
                TH(O) = THc*180/pi; PS(O) = Psi*180/pi;
            end
            if O == 4
                if std(TH)/mean(TH) < 0.001 && std(PS)/mean(PS) < 0.001; THc = mean(TH); Psi = mean(PS); ex = 1;
                else P = P + 1; TH = []; PS = []; O = 0;
                end
            end
        end
        TC = THc;
    elseif TYP == 20 || TYP == 30; TC = 0;
    elseif TYP == 25 || TYP == 35;
        Frang = SOLID.Frang*pi/180; SOLID.Flowr = OP.kc(S,T)*OP.R;
        eta = SOLID.mewL/tan(Frang); G = SOLID.Flowr/OP.R; ec = OP.R*(eta*(1-G)+(1-eta)*sqrt(1-G)); % ec - flow channel eccentricity (EN 1991-4:2006 Eqs 5.55-5.57)
        cosTHc = (OP.R*OP.R + ec*ec - SOLID.Flowr*SOLID.Flowr)/(2*OP.R*ec); THc = acos(cosTHc); % angular spread of wall contact with flow channel +- THc (EN 1991-4:2006 Eq. 5.58)
        Psi = asin((OP.R/SOLID.Flowr)*sin(THc));
        Uwc = 2*THc*OP.R; Usc = 2*SOLID.Flowr*(pi - Psi); % psi angle and arc lengths (EN 1991-4:2006 Eqs 5.59-5.61)
        Ac = (pi - Psi)*SOLID.Flowr*SOLID.Flowr + THc*OP.R*OP.R - OP.R*SOLID.Flowr*sin(Psi - THc); % Cross-sectional area of flow channel (EN 1991-4:2006 Eq. 5.62)
        zoc = (1/SOLID.K)*(Ac/(Uwc*SOLID.mewL + Usc*tan(Frang))); phco = SOLID.Weight*1e-6*SOLID.K*zoc; THc = THc*180/pi; % EN 1991-4:2006 Eqs 5.65-5.66
        OP.EN.THc = THc; OP.EN.Psi = Psi*180/pi; OP.EN.kc = OP.kc; OP.EN.ec = ec/OP.R; OP.EN.Ac = Ac/(pi*OP.R*OP.R); TC = THc;
    end
end

SET = 0; nS = 1; new = 0; SP = 0; SF = 0; DRAW3DS = 3;
for N = 1:length(NODES.index)
    if NODES.Z(N) <= OP.Hc % If not within roof
        if n == 0; n = 1; J = J + 1; Ntemp.index(J) = Noffset + J; Ntemp.R(J) = 0; Ntemp.T(J) = 0; Ntemp.Z(J) = 0;
        else
            if NODES.Z(N) ~= Ntemp.Z(J); new = 1;
                J = J + 1; Ntemp.index(J) = Noffset + J; Ntemp.R(J) = 0; Ntemp.T(J) = 0; Ntemp.Z(J) = NODES.Z(N);
            end
        end
        
        if SPRING.Discharge == 1
            if (NODES.T(N) >= TC) && (NODES.T(N) <= 360 - TC) && (OP.H - NODES.Z(N) > 0)
                S = S + 1; SPRING.index(S) = S; SPRING.nodes(S,1) = Ntemp.index(J); SPRING.nodes(S,2) = NODES.index(N); SPRING.region(S) = 2;

                tM = NODES.T(N); 
                if tM == 0; tR = NODES.T(N+1); tL = tR; 
                else 
                    if tM ~= OP.T; tR = NODES.T(N+1); else tR = NODES.T(N-1); end
                    tL = NODES.T(N-1); 
                end
                if tM == OP.T || tM == 360; tL = NODES.T(N-1); tR = tL; 
                else 
                    if tM ~= 0; tL = NODES.T(N-1); else tL = NODES.T(N+1); end
                    tR = NODES.T(N+1);
                end
                dTL = 0.5*abs(tL-tM); dTR = 0.5*abs(tM-tR);
                dT = dTL + dTR;

                zM = NODES.Z(N); zmid = OP.H - zM;
                if N + Ttot <= length(NODES.index); zU = NODES.Z(N+Ttot); else zU = NODES.Z(N-Ttot); end
                if N - Ttot >= 1; zB = NODES.Z(N-Ttot); else zB = NODES.Z(N+Ttot); end
                dZU = 0.5*abs(zU-zM); dZB = 0.5*abs(zM-zB);
                dZ = dZU + dZB;
                dA = dZ*dT*OP.R; % Shell area associated with each spring node

                % Calculation of k, spring stiffness
                if TYP == 19
                    Ac = (pi - Psi*pi/180)*SOLID.Flowr*SOLID.Flowr + pi*THc*OP.R*OP.R/180 - OP.R*SOLID.Flowr*sin((Psi - THc)*pi/180); As = pi*OP.R*OP.R - Ac;
                    Usw = 2*OP.R*(pi - pi*THc/180); Uwc = 2*THc*OP.R*pi/180; Usc = 2*SOLID.Flowr*(pi - Psi*pi/180);
                    zf = Ac/(Uwc*SOLID.K*SOLID.mewL + Usc*SOLID.K*tan(Frang)); % Eq. 9, within the channel, so lower mu (from EN)
                    zs = As/(Usw*SOLID.K*SOLID.mewU); % from Eq. 17, within the static solid, so upper mu (from EN)
                    u = zf/(zs + zf); n = tan(Frang)*SOLID.K*Usc*zf/As;
                    pV = SOLID.Weight*1e-6*zs*(1 + n + n*u*exp(-zmid/zf) - (1 + n + n*u)*exp(-zmid/zs)); % Eq. 17
                end
                Es = 3*(SOLID.Weight^(3/2))*pV; nCR = 4; vS = 0.3; b = 0.29*nCR;
                kk = (dA/OP.R)*b*Es*(1/(1-vS-2*vS*vS));

                if nS == 1
                    dA_old = dA; nS = 0; SET = SET + 1; SPRING.setlist.no(SET) = SET; SPRING.setlist.k(SET) = kk;
                else
                    if dA ~= dA_old || new == 1
                        SET = SET + 1; new = 0; SPRING.setlist.no(SET) = SET; SPRING.setlist.k(SET) = kk; end
                end
                SPRING.set(S) = SET; dA_old = dA; SPRING.dA(S) = dA;
                if DRAW3DS == 3
                    SP = SP + 1; ST_P.R(SP) = kk; ST_P.T(SP) = tM; ST_P.Z(SP) = zM;
                    SF = SF + 1; ST_F.R(SF) = dA; ST_F.T(SF) = tM; ST_F.Z(SF) = zM;
                end
            end
        else
            S = S + 1; SPRING.index(S) = S; SPRING.nodes(S,1) = Ntemp.index(J); SPRING.nodes(S,2) = NODES.index(N); SPRING.region(S) = 2;
        end
    end
end
if DRAW3DS == 3; AJSABQ_draw_3dsurf(ST_P,ST_F,OP,DRAW3DS,1); end

lenN = length(NODES.index);
lenT = length(Ntemp.index);
NODES.index(lenN+1:lenN+lenT) = Ntemp.index;
NODES.R(lenN+1:lenN+lenT) = Ntemp.R;
NODES.T(lenN+1:lenN+lenT) = Ntemp.T;
NODES.Z(lenN+1:lenN+lenT) = Ntemp.Z;
Axis = Ntemp.index;
clear('Ntemp');