function [fid,OP,ST_P,ST_F,SP,SF] = AJSABQ_fwrite_load(fid,S,T,OP,SOLID,ELEMENTS,NODES,QVV,DRAW3DS,DRAW2DP,ST_P,ST_F,SP,SF)
% Function to write the loading of the Step heding of the .inp file for the
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:24 (previously 22/10/11 - 16:42)

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

loading = OP.type(S,T); mag = OP.mag(S,T); factor = OP.factor(S,T); 
%ST_P = struct('R',[],'T',[],'Z',[]); ST_F = struct('R',[],'T',[],'Z',[]);
%SP = 0; SF = 0;
if loading == .1 || loading == .2 || loading == .3 || loading == .4 || loading == .5 || loading == .6 % Point load/moment through rigid point (bottom)
    disp('...writing point load/moment through rigid point (bottom)...');
    fprintf(fid,'%s\n','*CLOAD');
    fprintf(fid,'%s\n',['RIGID_NODE_BOTTOM, ',num2str(loading*10),', ',num2str(mag*factor)]);   
elseif loading == .7 || loading == .8 || loading == .9 % Point rotation through rigid point (bottom)
    disp('...writing point rotation through rigid point (bottom)...');
    fprintf(fid,'%s\n','*BOUNDARY, OP=NEW');
    fprintf(fid,'%s\n',['RIGID_NODE_BOTTOM, ',num2str(loading*10-3),', ',num2str(loading*10-3),', ',num2str(mag*factor)]);    
elseif loading == 1 || loading == 2 || loading == 3 || loading == 4 || loading == 5 || loading == 6 % Point load/moment through rigid point (top)
    disp('...writing point load/moment through rigid point (top)...');
    fprintf(fid,'%s\n','*CLOAD');
    fprintf(fid,'%s\n',['RIGID_NODE_TOP, ',num2str(loading),', ',num2str(mag*factor)]);   
elseif loading == 7 || loading == 8 || loading == 9 % Point rotation through rigid point (top)
    disp('...writing point rotation through rigid point (top)...');
    fprintf(fid,'%s\n','*BOUNDARY, OP=NEW');
    fprintf(fid,'%s\n',['RIGID_NODE_TOP, ',num2str(loading-3),', ',num2str(loading-3),', ',num2str(mag*factor)]);    
elseif loading == 10 % Cylinder top line load
    disp('...writing cylinder line load...');   
    if OP.Heltoggle == 1
        fprintf(fid,'%s\n','*DSLOAD, FOLLOWER=NO');
        fprintf(fid,'%s\n',['TOP_EDGE, EDLD, ',num2str(mag*factor),', 0., 0., -1.']);
    else
        fprintf(fid,'%s\n','*DSLOAD');
        fprintf(fid,'%s\n',['TOP_EDGE, EDNOR, ',num2str(mag*factor)]);
    end
elseif loading == 15 % Cylinder internal uniform pressure
    disp('...writing cylinder internal pressure...');
    fprintf(fid,'%s\n','*DSLOAD');
    fprintf(fid,'%s\n',['INNER_SURFACE_CYL, P, ',num2str(mag*factor)]);
elseif loading == 152 % Cylinder internal uniform friction
    disp('...writing cylinder internal friction...');
    fprintf(fid,'%s\n','*DSLOAD');
    fprintf(fid,'%s\n',['INNER_SURFACE_CYL, TRSHR, ',num2str(mag*factor),', 0., 0., -1.']);
elseif loading == 16 % Hopper internal uniform pressure
    disp('...writing hopper internal pressure...');
    fprintf(fid,'%s\n','*DSLOAD');
    fprintf(fid,'%s\n',['INNER_SURFACE_HOP, P, ',num2str(mag*factor)]);
elseif loading == 162 % Hopper internal uniform friction
    disp('...writing hopper internal friction...');
    fprintf(fid,'%s\n','*DSLOAD');
    fprintf(fid,'%s\n',['INNER_SURFACE_HOP, TRSHR, ',num2str(mag*factor),', 0., 0., -1.']);
elseif loading == 17 % Prescribed test circumferential/axial distribution
    disp('...writing prescribed inline function...');
    fprintf(fid,'%s\n',['** Definition of prescribed function: ',OP.Inline.S]);
    fprintf(fid,'%s\n','*DLOAD');
    for E = 1:length(ELEMENTS.index)-OP.Roof-OP.Hop
        e1 = ELEMENTS.nodes(E+OP.Hop,1); f1 = find(NODES.index == e1);
        e2 = ELEMENTS.nodes(E+OP.Hop,2); f2 = find(NODES.index == e2);
        e3 = ELEMENTS.nodes(E+OP.Hop,3); f3 = find(NODES.index == e3);
        Z1 = NODES.Z(f1); Z2 = NODES.Z(f2); Z3 = NODES.Z(f3);
        T1 = NODES.T(f1); T2 = NODES.T(f2); T3 = NODES.T(f3);
        zmid = OP.H - 0.5*(Z2 + Z3); Tcoord = 0.5*(T1 + T2);
        if T2 < T1 && Z3 > Z1; T2 = 360; T3 = 360; elseif Z3 < Z2 && T1 < T2; T1 = 360; end
        %if Tcoord >= OP.Inline.T0;
        if Tcoord <= OP.Inline.T0 || Tcoord >= 360 - OP.Inline.T0; %%%%% Change back to <= and >=
            pH = OP.Inline.F(Tcoord,OP.Inline.T0,zmid);
            fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E+OP.Hop),', P, ',num2str(pH*factor)]);
            if DRAW3DS == 1 || DRAW3DS == 3; SP = SP + 1; ST_P.R(SP) = pH; ST_P.T(SP) = Tcoord; ST_P.Z(SP) = zmid; end
        else
            if DRAW3DS == 1 || DRAW3DS == 3; SP = SP + 1; ST_P.R(SP) = 0; ST_P.T(SP) = Tcoord; ST_P.Z(SP) = zmid; end
        end
    end
    if DRAW3DS == 2 || DRAW3DS == 3; ST_F.R = -1; ST_F.T = -1; ST_F.Z = -1; end
elseif loading == 19 % Distributed Rotter 1986 eccentric silo pressures
    if loading == 19; disp(['...writing Rotter''s (1986) eccentric cylinder pressures (factored by ',num2str(factor),')...']); txt = ': ECCENTRIC'; wh = 'Rotter 1986'; end
    fprintf(fid,'%s\n',['** Definition of ',wh,' horizontal cylinder pressures and frictional tractions',txt,': (factored by ',num2str(factor),')']);
    fprintf(fid,'%s\n','*DLOAD'); 
    POS = struct('min2',[],'min',[],'max',[],'max2',[]);
    Frang = SOLID.Frang*pi/180; SOLID.Flowr = OP.kc(S,T)*OP.R;
    yk = tan(Frang)*SOLID.K/(SOLID.mewL*SOLID.K);
    vals = [OP.kc, yk]; tol = 0.01; N = 4; O = 0; P = 0; ex = 0; TH = []; PS = [];
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
        if O == 4;
            if std(TH)/mean(TH) < 0.001 && std(PS)/mean(PS) < 0.001; THc = mean(TH); Psi = mean(PS); ex = 1;
            else P = P + 1; TH = []; PS = []; O = 0;
            end
        end
    end
    Ac = (pi - Psi*pi/180)*SOLID.Flowr*SOLID.Flowr + pi*THc*OP.R*OP.R/180 - OP.R*SOLID.Flowr*sin((Psi - THc)*pi/180); As = pi*OP.R*OP.R - Ac;
    Usw = 2*OP.R*(pi - pi*THc/180); Uwc = 2*THc*OP.R*pi/180; Usc = 2*SOLID.Flowr*(pi - Psi*pi/180); 
    zf = Ac/(Uwc*SOLID.K*SOLID.mewL + Usc*SOLID.K*tan(Frang)); % Eq. 9, within the channel, so lower mu (from EN)
    pof = SOLID.K*SOLID.Weight*1e-6*zf; % Eq. 10
    zs = As/(Usw*SOLID.K*SOLID.mewU); % from Eq. 17, within the static solid, so upper mu (from EN)
    u = zf/(zs + zf); n = tan(Frang)*SOLID.K*Usc*zf/As;
    OP.Rot.THc = THc; OP.Rot.Psi = Psi; OP.Rot.kc = OP.kc; OP.Rot.ec = ecoR; OP.Rot.Ac = Ac/(pi*OP.R*OP.R);   
    POS.min2 = OP.chan0(S,T) - 2*THc; if POS.min2 < 0; POS.min2 = 360 + POS.min2; out = 1; end
    POS.min = OP.chan0(S,T) - THc; if POS.min < 0; POS.min = 360 + POS.min; out = 1; end
    POS.max = OP.chan0(S,T) + THc; if POS.max > 360; POS.max = POS.max - 360; out = 1; end
    POS.max2 = OP.chan0(S,T) + 2*THc; if POS.max2 > 360; POS.max2 = POS.max2 - 360; out = 1; end
    for E = 1:length(ELEMENTS.index)-OP.Roof-OP.Hop
        e1 = ELEMENTS.nodes(E+OP.Hop,1); f1 = find(NODES.index == e1);
        e2 = ELEMENTS.nodes(E+OP.Hop,2); f2 = find(NODES.index == e2);
        e3 = ELEMENTS.nodes(E+OP.Hop,3); f3 = find(NODES.index == e3);
        Z1 = NODES.Z(f1); Z2 = NODES.Z(f2); Z3 = NODES.Z(f3);
        T1 = NODES.T(f1); T2 = NODES.T(f2); T3 = NODES.T(f3);
        zmid = OP.H - 0.5*(Z2 + Z3); Tcoord = 0.5*(T1 + T2);
        if T2 < T1 && Z3 > Z1; T2 = 360; T3 = 360; elseif Z3 < Z2 && T1 < T2; T1 = 360; end
        % Static pressures
        ph = SOLID.K*SOLID.Weight*1e-6*zs*(1 + n + n*u*exp(-zmid/zf) - (1 + n + n*u)*exp(-zmid/zs)); % Eq. 17     
        pf = SOLID.mewU*ph; pff = pf; phh = ph; 
        % Flow channel pressures
        phce = pof*(1 - exp(-zmid/zf)); % Eq. 8 or 11
        phae = 2*ph - phce; pfce = SOLID.mewL*phce; pfae = SOLID.mewL*(2*pf - pfce);
        if E == 1; OP.Rot.phs = ph*1e3; OP.Rot.phc = phce*1e3; end
        if Tcoord >= POS.min2 && Tcoord < POS.min % within the right adjacent region
            if OP.ears(S,T) == 1; phh = phae; pff = pfae; else phh = ph; pff = pf; end
        elseif (out == 0 && Tcoord >= POS.min && Tcoord < POS.max) % within the flow channel
            phh = phce; pff = pfce;
        elseif ((out == 1 && Tcoord >= POS.min && Tcoord <= 360) || (out == 1 && Tcoord >= 0 && Tcoord < POS.max));
            phh = phce; pff = pfce;
        elseif Tcoord >= POS.max && Tcoord < POS.max2 % within the left adjacent region
            if OP.ears(S,T) == 1; phh = phae; pff = pfae; else phh = ph; pff = pf; end
        else phh = ph; pff = pf;
        end
        fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E+OP.Hop),', P, ',num2str(phh*factor)]);
        fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E+OP.Hop),', TRSHR, ',num2str(pff*factor),', 0., 0., -1.']);
        if DRAW3DS == 1 || DRAW3DS == 3; SP = SP + 1; ST_P.R(SP) = phh*factor; ST_P.T(SP) = Tcoord; ST_P.Z(SP) = 0.5*(Z2+Z3); end
        if DRAW3DS == 2 || DRAW3DS == 3; SF = SF + 1; ST_F.R(SF) = pff*factor; ST_F.T(SF) = Tcoord; ST_F.Z(SF) = 0.5*(Z2+Z3); end
    end
    if DRAW2DP == 1; AJSABQ_draw_2dplanar(OP.Rot,S); end
elseif loading == 20 || loading == 201 || loading == 25 || loading == 30 || loading == 35 % Distributed Janssen/modified Reimbert concentric/eccentric silo pressures
     if loading == 20 || loading == 201 || loading == 30 % Apply EN 1991-4 discharge factors or not for concentric discharge only
        if OP.EN_Dis(S,T) == 1
            dis_Ch = SOLID.Ch; dis_Cw = SOLID.Cw; 
            txEN = [': EN 1991-4 discharge factors Ch = ',num2str(SOLID.Ch),'; Cw = ',num2str(SOLID.Cw)]; txE = ' and EN 1991-4 Ch & Cw';
        else   
            dis_Ch = 1; dis_Cw = 1; txEN = ''; txE = '';
        end
    else    
        dis_Ch = 1; dis_Cw = 1; txEN = ''; txE = '';
    end    
    if loading == 20; disp(['...writing Janssen concentric cylinder pressures (factored by ',num2str(factor),txE,')...']); txt = ': CONCENTRIC'; wh = 'Janssen'; end
    if loading == 201; disp(['...writing Janssen concentric cylinder pressures + custom inline load (factored by ',num2str(factor),txE,')...']); txt = ': CONCENTRIC'; wh = 'Janssen'; end
    if loading == 25; disp(['...writing Janssen eccentric cylinder pressures (factored by ',num2str(factor),')...']); txt = ': ECCENTRIC'; wh = 'Janssen'; end
    if loading == 30; disp(['...writing mod. Reimbert concentric cylinder pressures (factored by ',num2str(factor),txE,')...']); txt = ': CONCENTRIC'; wh = 'modified Reimbert'; end
    if loading == 35; disp(['...writing mod. Reimbert eccentric cylinder pressures (factored by ',num2str(factor),')...']); txt = ': ECCENTRIC'; wh = 'modified Reimbert'; end
    fprintf(fid,'%s\n',['** Definition of ',wh,' horizontal cylinder pressures and frictional tractions',txt,': (factored by ',num2str(factor),')',txEN]);
    fprintf(fid,'%s\n','*DLOAD'); out = 0;
    if loading == 25 || loading == 35
        POS = struct('min2',[],'min',[],'max',[],'max2',[]);
        Frang = SOLID.Frang*pi/180; SOLID.Flowr = OP.kc(S,T)*OP.R;
        eta = SOLID.mewL/tan(Frang); G = SOLID.Flowr/OP.R; ec = OP.R*(eta*(1-G)+(1-eta)*sqrt(1-G)); % ec - flow channel eccentricity (EN 1991-4:2006 Eqs 5.55-5.57)
        cosTHc = (OP.R*OP.R + ec*ec - SOLID.Flowr*SOLID.Flowr)/(2*OP.R*ec); THc = acos(cosTHc); % angular spread of wall contact with flow channel +- THc (EN 1991-4:2006 Eq. 5.58)
        Psi = asin((OP.R/SOLID.Flowr)*sin(THc));
        Uwc = 2*THc*OP.R; Usc = 2*SOLID.Flowr*(pi - Psi); % psi angle and arc lengths (EN 1991-4:2006 Eqs 5.59-5.61)
        Ac = (pi - Psi)*SOLID.Flowr*SOLID.Flowr + THc*OP.R*OP.R - OP.R*SOLID.Flowr*sin(Psi - THc); % Cross-sectional area of flow channel (EN 1991-4:2006 Eq. 5.62)
        zoc = (1/SOLID.K)*(Ac/(Uwc*SOLID.mewL + Usc*tan(Frang))); phco = SOLID.Weight*1e-6*SOLID.K*zoc; THc = THc*180/pi; % EN 1991-4:2006 Eqs 5.65-5.66
        OP.EN.THc = THc; OP.EN.Psi = Psi*180/pi; OP.EN.kc = OP.kc; OP.EN.ec = ec/OP.R; OP.EN.Ac = Ac/(pi*OP.R*OP.R);
        POS.min2 = OP.chan0(S,T) - 2*THc; if POS.min2 < 0; POS.min2 = 360 + POS.min2; out = 1; end
        POS.min = OP.chan0(S,T) - THc; if POS.min < 0; POS.min = 360 + POS.min; out = 1; end
        POS.max = OP.chan0(S,T) + THc; if POS.max > 360; POS.max = POS.max - 360; out = 1; end
        POS.max2 = OP.chan0(S,T) + 2*THc; if POS.max2 > 360; POS.max2 = POS.max2 - 360; out = 1; end
    end
    for E = 1:length(ELEMENTS.index)-OP.Roof-OP.Hop
        e1 = ELEMENTS.nodes(E+OP.Hop,1); f1 = find(NODES.index == e1);
        e2 = ELEMENTS.nodes(E+OP.Hop,2); f2 = find(NODES.index == e2);
        e3 = ELEMENTS.nodes(E+OP.Hop,3); f3 = find(NODES.index == e3);
        Z1 = NODES.Z(f1); Z2 = NODES.Z(f2); Z3 = NODES.Z(f3);
        T1 = NODES.T(f1); T2 = NODES.T(f2); T3 = NODES.T(f3);
        zmid = OP.H - 0.5*(Z2 + Z3); Tcoord = 0.5*(T1 + T2);
        if T2 < T1 && Z3 > Z1; T2 = 360; T3 = 360; elseif Z3 < Z2 && T1 < T2; T1 = 360; end
        % Static pressures
        if loading == 20 || loading == 201 || loading == 25; ph = OP.Con.pho*(1 - exp(-zmid/OP.Con.zo));
        elseif loading == 30 || loading == 35; ph = OP.Con.pho*(1 - ((zmid - SOLID.Equiv)/(OP.Con.zo - SOLID.Equiv) + 1)^OP.Con.nR);
        end
        if loading == 201
            if Tcoord <= OP.Inline.T0 || Tcoord >= 360 - OP.Inline.T0
                ph = ph + ph*OP.Inline.F(Tcoord,OP.Inline.T0,zmid); 
            end
        end
        pf = SOLID.mewU*ph; pff = pf*dis_Cw; phh = ph*dis_Ch;
        if loading == 25 || loading == 35
            % Flow channel pressures
            phce = phco*(1 - exp(-zmid/zoc)); pfce = SOLID.mewL*phce;
            phae = 2*ph - phce;  pfae = SOLID.mewL*phae; 
            if Tcoord >= POS.min2 && Tcoord < POS.min % within the right adjacent region
                if OP.ears(S,T) == 1; phh = phae; pff = pfae; else phh = ph; pff = pf; end
            elseif (out == 0 && Tcoord >= POS.min && Tcoord < POS.max) % within the flow channel
                phh = phce; pff = pfce;
            elseif ((out == 1 && Tcoord >= POS.min && Tcoord <= 360) || (out == 1 && Tcoord >= 0 && Tcoord < POS.max))
                phh = phce; pff = pfce;
            elseif Tcoord >= POS.max && Tcoord < POS.max2 % within the left adjacent region
                if OP.ears(S,T) == 1; phh = phae; pff = pfae; else phh = ph; pff = pf; end
            else phh = ph; pff = pf;
            end
            if E == 1; OP.EN.phs = ph*1e3; OP.EN.phc = phce*1e3; end
        end
        fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E+OP.Hop),', P, ',num2str(phh*factor)]); 
        fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E+OP.Hop),', TRSHR, ',num2str(pff*factor),', 0., 0., -1.']); 
        if DRAW3DS == 1 || DRAW3DS == 3; SP = SP + 1; ST_P.R(SP) = phh*factor; ST_P.T(SP) = Tcoord; ST_P.Z(SP) = 0.5*(Z2+Z3); end
        if DRAW3DS == 2 || DRAW3DS == 3; SF = SF + 1; ST_F.R(SF) = pff*factor; ST_F.T(SF) = Tcoord; ST_F.Z(SF) = 0.5*(Z2+Z3); end       
    end
    if loading == 25 && DRAW2DP == 1 || loading == 35 && DRAW2DP == 1; AJSABQ_draw_2dplanar(OP.EN,S); end
elseif loading == 40 || loading == 45 || loading == 46 % Distributed Walker concentric hopper pressures, plus eccentric depression
    if loading == 40; disp(['...writing Walker concentric hopper pressures (factored by ',num2str(factor),')...']); txt = ': CONCENTRIC'; end
    if loading == 45; disp(['...writing Walker eccentric C1 hopper pressures (factored by ',num2str(factor),')...']); txt = ': C1 ECCENTRIC'; end
    if loading == 46; disp(['...writing Walker eccentric C2 hopper pressures (factored by ',num2str(factor),')...']); txt = ': C2 ECCENTRIC'; end
    fprintf(fid,'%s\n',['** Definition of Walker normal hopper pressures and frictional tractions',txt,': (factored by ',num2str(factor),')']);
    fprintf(fid,'%s\n','*DLOAD'); th0 = OP.hop_th0(S,T)*pi/180; eta = OP.hop_eta(S,T);
    for E = 1:OP.Hop
        e1 = ELEMENTS.nodes(E,1); f1 = find(NODES.index == e1);
        e2 = ELEMENTS.nodes(E,2); f2 = find(NODES.index == e2);
        e3 = ELEMENTS.nodes(E,3); f3 = find(NODES.index == e3);
        Z1 = NODES.Z(f1); Z2 = NODES.Z(f2); Z3 = NODES.Z(f3);
        T1 = NODES.T(f1); T2 = NODES.T(f2); T3 = NODES.T(f3);
        R1 = NODES.R(f1); R2 = NODES.R(f2); R3 = NODES.R(f3);
        zmid = 0.5*(Z2 + Z3); Tc = 0.5*(T2 + T1)*pi/180;
        if T2 < T1 && Z3 > Z1; T2 = 360; T3 = 360;
        elseif Z3 < Z2 && T1 < T2; T1 = 360;
        end
        pv = (SOLID.Weight*1e-6*OP.Hh)/(OP.ConH.n - 1)*((zmid/OP.Hh) - (zmid/OP.Hh)^OP.ConH.n) + OP.ConH.pvft*(zmid/OP.Hh)^OP.ConH.n;
        if loading == 45 % C1 patch ('smooth')
            if Tc <= th0 || Tc >= (2*pi - th0); pv = pv + 0.5*(1 + cos(pi*Tc/th0))*eta*pv; end
        elseif loading == 46 % C2 patch ('hypersmooth')
            if TC <= th0 || Tc >= (2*pi - th0); pv = pv + (cos(0.5*pi*Tc/th0))*(cos(0.5*pi*Tc/th0))*(cos(0.5*pi*Tc/th0))*(cos(0.5*pi*Tc/th0))*eta*pv; end
        end
        pn = OP.ConH.F*pv; pf = pn*SOLID.mewU;
        fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E),', P, ',num2str(pn*factor)]);
        fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E),', TRSHR,' ,num2str(pf*factor),', 0., 0., -1.']);
        if DRAW3DS == 1 || DRAW3DS == 3; SP = SP + 1; ST_P.R(SP) = pn*factor; ST_P.T(SP) = Tcoord; ST_P.Z(SP) = zmid; end
        if DRAW3DS == 2 || DRAW3DS == 3; SF = SF + 1; ST_F.R(SF) = pf*factor; ST_F.T(SF) = Tcoord; ST_F.Z(SF) = zmid; end
    end
elseif loading == 50 % EN 1993-4-1 Appendix C Wind pressure distribution (single silo)
    disp('...writing EN 1993-4-1 Appendix C cylinder wind pressures (single silo)...');
    fprintf(fid,'%s\n',['** Definition of EN 1993-4-1 Appendix C cylinder wind pressures (single silo) with stagnation pressure ',num2str(mag),', (factored by ',num2str(factor),')']);
    fprintf(fid,'%s\n','*DLOAD');
    for E = 1:length(ELEMENTS.index)-OP.Roof-OP.Hop
        e1 = ELEMENTS.nodes(E+OP.Hop,1); f1 = find(NODES.index == e1);
        e2 = ELEMENTS.nodes(E+OP.Hop,2); f2 = find(NODES.index == e2);
        e3 = ELEMENTS.nodes(E+OP.Hop,3); f3 = find(NODES.index == e3);
        Z1 = NODES.Z(f1); Z2 = NODES.Z(f2); Z3 = NODES.Z(f3);
        T1 = NODES.T(f1); T2 = NODES.T(f2); T3 = NODES.T(f3);
        zmid = OP.H - 0.5*(Z2 + Z3); Tcoord = 0.5*(T1 + T2);
        if T2 < T1 && Z3 > Z1; T2 = 360; T3 = 360; elseif Z3 < Z2 && T1 < T2; T1 = 360; end
        doH = 2*OP.R/OP.Hc; % Aspect ratio for complete structure (should strictly include supports!)
        Cp = -0.54 + 0.16*doH + (0.28 + 0.04*doH)*cos(Tcoord*pi/180) + (1.04 - 0.20*doH)*cos(2*Tcoord*pi/180)...
            + (0.36 - 0.05*doH)*cos(3*Tcoord*pi/180) - (0.14 - 0.05*doH)*cos(4*Tcoord*pi/180); % Eq. C.1 - assuming closed roof
        fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E),', P, ',num2str(Cp*mag*factor)]);
        if DRAW3DS == 1 || DRAW3DS == 3; SP = SP + 1; ST_P.R(SP) = Cp*mag*factor; ST_P.T(SP) = Tcoord; ST_P.Z(SP) = zmid; end
    end  
elseif loading == 51 % EN 1993-4-1 Appendix C Wind pressure distribution (group silo)
    disp('...writing EN 1993-4-1 Appendix C cylinder wind pressures (group silo)...');
    fprintf(fid,'%s\n',['** Definition of EN 1993-4-1 Appendix C cylinder wind pressures *group silo) with stagnation pressure ',num2str(mag),', (factored by ',num2str(factor),')']);
    fprintf(fid,'%s\n','*DLOAD');
    for E = 1:length(ELEMENTS.index)-OP.Roof-OP.Hop;
        e1 = ELEMENTS.nodes(E+OP.Hop,1); f1 = find(NODES.index == e1);
        e2 = ELEMENTS.nodes(E+OP.Hop,2); f2 = find(NODES.index == e2);
        e3 = ELEMENTS.nodes(E+OP.Hop,3); f3 = find(NODES.index == e3);
        T1 = NODES.T(f1); T2 = NODES.T(f2); T3 = NODES.T(f3); Tcoord = 0.5*(T1 + T2);
        if T2 < T1 && Z3 > Z1; T2 = 360; T3 = 360; elseif Z3 < Z2 && T1 < T2; T1 = 360; end
        doH = 2*OP.R/OP.Hc; % Aspect ratio for complete structure (should strictly include supports!)
        Cp = 0.20 + 0.60*cos(Tcoord*pi/180) + 0.27*cos(2*Tcoord*pi/180) - 0.05*cos(3*Tcoord*pi/180) - 0.13*cos(4*Tcoord*pi/180)...
            + 0.13*cos(6*Tcoord*pi/180) - 0.09*cos(8*Tcoord*pi/180) + 0.07*cos(10*Tcoord*pi/180); % Eq. C.2 - assuming closed roof
        fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E),', P, ',num2str(Cp*mag*factor)]);
        if DRAW3DS == 1 || DRAW3DS == 3; SP = SP + 1; ST_P.R(SP) = Cp*mag*factor; ST_P.T(SP) = Tcoord; ST_P.Z(SP) = zmid; end
    end
elseif loading == 100 % 100 = External silo pressure data from .cvv file
    disp('...writing externally generated .cvv file cylinder pressures...');
    fprintf(fid,'%s\n',['** Definition of external silo pressures and friction tractions from .cvv file: ',OP.cvv_file]);
    fprintf(fid,'%s\n','*DLOAD');
    for E = 1:length(ELEMENTS.index)-OP.Roof-OP.Hop
        e1 = ELEMENTS.nodes(E+OP.Hop,1); f1 = find(NODES.index == e1);
        e2 = ELEMENTS.nodes(E+OP.Hop,2); f2 = find(NODES.index == e2);
        e3 = ELEMENTS.nodes(E+OP.Hop,3); f3 = find(NODES.index == e3);
        Z1 = NODES.Z(f1)-OP.Hh; Z2 = NODES.Z(f2)-OP.Hh; Z3 = NODES.Z(f3)-OP.Hh;
        T1 = NODES.T(f1); T2 = NODES.T(f2); T3 = NODES.T(f3);
        if T2 < T1 && Z3 > Z1; T2 = 360; T3 = 360; elseif Z3 < Z2 && T1 < T2; T1 = 360; end
        t_mid = 0.5*(T1 + T2); z_mid = 0.5*(Z2 + Z3); dTt = (OP.T/QVV.dTT+1); goZ = 1;
        for TZ = 1:dTt:length(QVV.Z)-dTt
            z1 = QVV.Z(TZ); z2 = QVV.Z(TZ+dTt);
            if (z_mid >= z1 && z_mid <= z2 && goZ == 1) || (z_mid <= z1 && z_mid >= z2 && goZ == 1)
                z_frac = abs(z_mid-z1)/abs(z2-z1); goZ = 0; goT = 1;
                for TT = 0:1:dTt-1
                    t1 = QVV.T(TZ+TT); t2 = QVV.T(TZ+TT+1);
                    if (t_mid >= t1 && t_mid <= t2 && goT == 1)
                        t_frac = abs(t_mid-t1)/abs(t2-t1); goT = 0;
                        q1 = QVV.qn(TZ+TT); q2 = QVV.qn(TZ+TT+1);
                        q3 = QVV.qn(TZ+dTt+TT); q4 = QVV.qn(TZ+dTt+TT+1);
                        qn_tbot = q1 + t_frac*(q2-q1); qn_ttop = q3 + t_frac*(q4-q3); qn_t = 0.5*(qn_tbot+qn_ttop);
                        qn_zleft = q1 + z_frac*(q3-q1); qn_zright = q2 + z_frac*(q4-q2); qn_z = 0.5*(qn_zleft+qn_zright);
                        phh = 0.5*(qn_z+qn_t); pff = phh*QVV.SOLID.mu_w;
                        fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E+OP.Hop),', P, ',num2str(phh*factor)]);
                        fprintf(fid,'%s\n',['SHELL_INSTANCE.',num2str(E+OP.Hop),', TRSHR, ',num2str(pff*factor),', 0., 0., -1.']);
                        if DRAW3DS == 1 || DRAW3DS == 3; SP = SP + 1; ST_P.R(SP) = phh*factor; ST_P.T(SP) = t_mid; ST_P.Z(SP) = z_mid; end
                        if DRAW3DS == 2 || DRAW3DS == 3; SF = SF + 1; ST_F.R(SF) = pff*factor; ST_F.T(SF) = t_mid; ST_F.Z(SF) = z_mid; end
                    end
                end
            end
        end
    end
end