function [GEO,RSL,STEP,OP,NAME] = AJSABQ_assemble(GEO,RSL,STEP,SOLID,OP,NAME,MATERIAL)
% Function to check the assemble the resolution struct for the 
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 12:47 (previously 23/02/12 - 15:24)

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

OP.R = GEO.Rtot; OP.Sol = 0; GEO.Ttot = RSL.CIRC.t2(length(RSL.CIRC.t2)); OP.T = GEO.Ttot; I = 0;

if GEO.hop.Inc ~= 1; RSL.HOP.h1 = []; RSL.HOP.h2 = []; RSL.HOP.Z = []; RSL.HOP.t = []; RSL.HOP.Type = []; end
if GEO.cyl.Inc ~= 1; RSL.CYL.h1 = []; RSL.CYL.h2 = []; RSL.CYL.Z = []; RSL.CYL.t = []; RSL.CYL.Type = []; end

if ~isempty(RSL.HOP.h1) && GEO.hop.Inc == 1
    for H = 1:length(RSL.HOP.h1); I = I + 1;
        RSL.Region(I) = I; RSL.Where(I) = 1; RSL.h1(I) = RSL.HOP.h1(H); RSL.h2(I) = RSL.HOP.h2(H); 
        RSL.Z(I) = RSL.HOP.Z(H); RSL.t(I) = RSL.HOP.t(H); RSL.Type(I) = RSL.HOP.Type(H);
    end
end

if ~isempty(RSL.CYL.h1) && GEO.cyl.Inc == 1
    for H = 1:length(RSL.CYL.h1); I = I + 1;
        RSL.Region(I) = I; RSL.Where(I) = 2; RSL.h1(I) = RSL.CYL.h1(H); RSL.h2(I) = RSL.CYL.h2(H); 
        RSL.Z(I) = RSL.CYL.Z(H); RSL.t(I) = RSL.CYL.t(H); RSL.Type(I) = RSL.CYL.Type(H);
    end
end


RSL.t1 = RSL.CIRC.t1; RSL.t2 = RSL.CIRC.t2; RSL.T = RSL.CIRC.T; RSL.RoofZ = RSL.CYL.RoofZ; mem = []; GEO.Ttot = RSL.CIRC.t2(length(RSL.CIRC.t2));
for I = 1:length(RSL.Where)
    if GEO.hop.Inc ~= 1 && RSL.Where(I) == 1; mem(length(mem)+1) = I; end
    if GEO.cyl.Inc ~= 1 && RSL.Where(I) == 2; mem(length(mem)+1) = I; end
end
RSL.Where(mem) = []; RSL.h1(mem) = []; RSL.h2(mem) = []; RSL.t(mem) = []; RSL.Type(mem) = []; RSL.Z(mem) = []; RSL.Ttot = sum(RSL.T);
if GEO.cyl.Inc == 1; OP.Hc = GEO.cyl.Htot; else OP.Hc = 0; end
if GEO.hop.Inc == 1; OP.Hh = GEO.hop.Htot; else OP.Hh = 0; end
OP.H = OP.Hc + OP.Hh; OP.Htot = OP.H + tan(pi*GEO.roof.Angle/180)*OP.R;

if GEO.roof.Inc == 1; L = length(RSL.Where)+1; % Roof data
    RSL.Region(L) = L; RSL.Where(L) = 3; RSL.h1(L) = OP.Hc + OP.Hh; RSL.h2(L) = OP.Hc + OP.Hh + tan(pi*GEO.roof.Angle/180)*OP.R;
    RSL.Z(L) = RSL.RoofZ; RSL.t(L) = GEO.roof.Thick; RSL.Type(L) = 'c';
end
RSL.Type = char(RSL.Type);

if GEO.feature.Super_ell(1) ~= 1; GEO.feature.Super_ell = 0; end
if GEO.feature.Axi_sin(1) ~= 1; GEO.feature.Axi_sin = 0; end
if GEO.feature.Circ_cos(1) ~= 1; GEO.feature.Circ_cos = 0; end

% NOTE - Element types; a - S4; b - S4R; c - S4R5; d - S3; e - STRI3; f - STRI65; g - S8R; h - S8R5; j - S9R5
RSL.Hoptot = 0; RSL.Cyltot = 0; RSL.Rooftot = 0;
for I = 1:length(RSL.Where)
    if RSL.Type(I) == 'a' || RSL.Type(I) == 'b' || RSL.Type(I) == 'c' % 4-node elements
        RSL.Els(I) = RSL.Z(I)*RSL.Ttot;
    elseif RSL.Type(I) == 'd' || RSL.Type(I) == 'e' % 3-node elements
        RSL.Els(I) = RSL.Z(I)*RSL.Ttot*2;
    elseif RSL.Type(I) == 'f' % 6-node element
        RSL.Els(I) = RSL.Z(I)*RSL.Ttot/2;
    elseif RSL.Type(I) == 'g' || RSL.Type(I) == 'h' || RSL.Type(I) == 'j' % 8-node and 9-node elements
        RSL.Els(I) = RSL.Z(I)*RSL.Ttot/4;
    end
    if GEO.hop.Inc == 1 && RSL.Where(I) == 2 % In cylinder
        RSL.h1(I) = RSL.h1(I) + GEO.hop.Htot; RSL.h2(I) = RSL.h2(I) + GEO.hop.Htot;
    end
    if RSL.Where(I) == 1; RSL.Hoptot = RSL.Hoptot + RSL.Els(I); end
    if RSL.Where(I) == 2; RSL.Cyltot = RSL.Cyltot + RSL.Els(I); end
    if RSL.Where(I) == 3; RSL.Rooftot = RSL.Rooftot + RSL.Els(I); end
end

OP.Hop = RSL.Hoptot; OP.Cyl = RSL.Cyltot; OP.Roof = RSL.Rooftot;

%RSL.Region = [1:length(RSL.h1)]; 
if GEO.hop.Inc == 1; GEO.hop.Beta = atan(GEO.Rtot/GEO.hop.Htot)*180/pi; else GEO.hop.Beta = 0; end
RSL.Etot = sum(RSL.Els); STEP.Tot = length(STEP.Ana.Type);

I = 0; J = 0; S = 0; K = 0;
temp = STEP.Load.Ears; STEP.Load.Ears = []; 
temp2 = STEP.Load.Factor; STEP.Load.Factor = [];
temp3 = STEP.Load.EN_Dis; STEP.Load.EN_Dis = []; 
for L = 1:size(STEP.Load.Type,1)
    % Ears are ON by default, the following sets them to off wherever it was requested
    if STEP.Load.Type(L,3) == 19 || STEP.Load.Type(L,3) == 20 || STEP.Load.Type(L,3) == 201 || STEP.Load.Type(L,3) == 25 || STEP.Load.Type(L,3) == 30 || STEP.Load.Type(L,3) == 35; 
        I = I + 1;
        STEP.Load.Ears(I,1) = STEP.Load.Type(L,1); 
        STEP.Load.Ears(I,2) = STEP.Load.Type(L,2); 
        STEP.Load.Ears(I,3) = 1;
        for X = 1:size(temp,1)
            if STEP.Load.Ears(I,1) == temp(X,1) && STEP.Load.Ears(I,2) == temp(X,2); STEP.Load.Ears(I,3) = 0; end
        end
    end
    
    % EN 1991-4 discharge factors are OFF by default, the following sets them to on wherever it was requested
    K = K + 1;
    STEP.Load.EN_Dis(K,1) = STEP.Load.Type(L,1);
    STEP.Load.EN_Dis(K,2) = STEP.Load.Type(L,2);
    STEP.Load.EN_Dis(K,3) = 0;
    if STEP.Load.Type(L,3) == 20 || STEP.Load.Type(L,3) == 201 || STEP.Load.Type(L,3) == 30
        for X = 1:size(temp3,1)
            if STEP.Load.EN_Dis(K,1) == temp3(X,1) && STEP.Load.EN_Dis(K,2) == temp3(X,2); STEP.Load.EN_Dis(K,3) = 1; end
        end   
    end    
    
    STEP.Load.Factor(L,1) = STEP.Load.Type(L,1); STEP.Load.Factor(L,2) = STEP.Load.Type(L,2); STEP.Load.Factor(L,3) = 1;
    for X = 1:size(temp2,1)
        if STEP.Load.Factor(L,1) == temp2(X,1) && STEP.Load.Factor(L,2) == temp2(X,2); STEP.Load.Factor(L,3) = temp2(X,3); end
    end
    if L == 1 && size(STEP.Load.Type,1) == 1; STEP.Sub(1) = 1; end
    if L == 1; J = J + 1; continue
    else
        if STEP.Load.Type(L,1) == STEP.Load.Type(L-1,1)
            if L < size(STEP.Load.Type,1); J = J + 1; end
            if L == size(STEP.Load.Type,1); J = J + 1; S = S + 1; STEP.Sub(S) = J; end
        elseif STEP.Load.Type(L,1) ~= STEP.Load.Type(L-1,1)
            if L < size(STEP.Load.Type,1); S = S + 1; STEP.Sub(S) = J; J = 1; end
            if L == size(STEP.Load.Type,1); S = S + 1; STEP.Sub(S) = J; S = S + 1; STEP.Sub(S) = 1; end
        end
    end
end

if GEO.hop.Inc == 1; temp = STEP.BC.Pit; STEP.BC.Pit = []; 
    for L = 1:length(STEP.Ana.Type)
        STEP.BC.Pit(L,1) = L; STEP.BC.Pit(L,2:8) = zeros(1,7);        
        for X = 1:size(temp,1)
            if STEP.BC.Pit(L,1) == temp(X,1); STEP.BC.Pit(L,2:8) = temp(X,2:8); end
        end
    end
end

if GEO.cyl.Inc == 1
    temp = STEP.BC.Top; STEP.BC.Top = [];
    for L = 1:length(STEP.Ana.Type)
        STEP.BC.Top(L,1) = L; STEP.BC.Top(L,2:7) = zeros(1,6);        
        for X = 1:size(temp,1)
            if STEP.BC.Top(L,1) == temp(X,1); STEP.BC.Top(L,2:7) = temp(X,2:7); end
        end
    end
    temp = STEP.BC.Top_RigidRing; STEP.BC.Top_RigidRing = [];
    for L = 1:length(STEP.Ana.Type)
        STEP.BC.Top_RigidRing(L,1) = L; STEP.BC.Top_RigidRing(L,2:7) = zeros(1,6);        
        for X = 1:size(temp,1)
            if STEP.BC.Top_RigidRing(L,1) == temp(X,1); STEP.BC.Top_RigidRing(L,2:7) = temp(X,2:7); end
        end
    end   
end

temp = STEP.BC.Trans; STEP.BC.Trans = [];
for L = 1:length(STEP.Ana.Type)
    STEP.BC.Trans(L,1) = L; STEP.BC.Trans(L,2:7) = zeros(1,6);
    for X = 1:size(temp,1)
        if STEP.BC.Trans(L,1) == temp(X,1); STEP.BC.Trans(L,2:7) = temp(X,2:7); end
    end
end
temp = STEP.BC.Trans_RigidRing; STEP.BC.Trans_RigidRing = [];
for L = 1:length(STEP.Ana.Type)
    STEP.BC.Trans_RigidRing(L,1) = L; STEP.BC.Trans_RigidRing(L,2:7) = zeros(1,6);
    for X = 1:size(temp,1)
        if STEP.BC.Trans_RigidRing(L,1) == temp(X,1); STEP.BC.Trans_RigidRing(L,2:7) = temp(X,2:7); end
    end
end

if OP.T == 360
    STEP.BC.T0(:,2) = zeros(1,STEP.Tot); 
    STEP.BC.Ttot(:,2) = zeros(1,STEP.Tot); 
end

if ~isempty(find(STEP.Load.Type(:,3) == 100, 1)) 
    OP.cvv = 1; OP.cvv_file = STEP.Load.QVV; OP.cvv_file(length(OP.cvv_file)+1:length(OP.cvv_file)+4) = '.qvv';
else 
    OP.cvv = 0; 
end

OP.CWeld = 0;
if GEO.feature.C_Weld.Toggle ~= 1; GEO.feature.C_Weld.Type = []; GEO.feature.C_Weld.Z = []; GEO.feature.C_Weld.Ampt = []; end
if GEO.feature.C_Weld.Toggle == 1; GEO.feature.C_Weld.Z = GEO.feature.C_Weld.Z + OP.Hh; OP.CWeld = 1; end

if GEO.feature.A_Weld.Toggle ~= 1; GEO.feature.A_Weld.Z = []; GEO.feature.A_Weld.Ampt = []; end
if GEO.feature.A_Weld.Toggle == 1; GEO.feature.A_Weld.Z = GEO.feature.A_Weld.Z + OP.Hh; end

if isempty(GEO.feature.LBA) == 1; GEO.feature.LBA = 0; end
if isempty(GEO.feature.Circ_cos) == 1; GEO.feature.Circ_cos = 0; end
if isempty(GEO.feature.Axi_sin) == 1; GEO.feature.Axi_sin = 0; end
if isempty(GEO.feature.Plate_level) == 1; GEO.feature.Plate_level = 0; end

if ~isempty(find(STEP.Load.Type(:,3) == 20, 1)) || ~isempty(find(STEP.Load.Type(:,3) == 201, 1)) ...
        || ~isempty(find(STEP.Load.Type(:,3) == 25, 1)) || ~isempty(find(STEP.Load.Type(:,3) == 30, 1)) ...
        || ~isempty(find(STEP.Load.Type(:,3) == 35, 1))
    % Janssen/modified Reimbert concentric static soil pressures (EN 1991-4:2006 Eqs 5.1-5.6 or 5.71 - 5.76)
    OP.Sol = 1; OP.Con.Ar = pi*OP.R*OP.R; OP.Con.U = 2*pi*OP.R; OP.Con.zo = OP.Con.Ar/(SOLID.K*SOLID.mewU*OP.Con.U);
    OP.Con.pho = SOLID.Weight*1e-6*SOLID.K*OP.Con.zo; OP.Con.nR = -(1 + tan(SOLID.Repose*pi/180))*(1-SOLID.Equiv/OP.Con.zo); 
end
if ~isempty(find(STEP.Load.Type(:,3) == 40, 1)) || ~isempty(find(STEP.Load.Type(:,3) == 45, 1)) || ~isempty(find(STEP.Load.Type(:,3) == 46, 1))
    % Walker/Dabrowski hopper pressures
    OP.ConH.F = (1 + 0.8*SOLID.mewU*cot(GEO.hop.Beta*pi/180)); OP.ConH.n = 2*(OP.ConH.F*SOLID.mewU*cot(GEO.hop.Beta*pi/180) + OP.ConH.F - 1);
    OP.ConH.pvft = SOLID.Weight*1e-6*OP.R;   
end

NAM = NAME.FileName;
NAME.FileName(length(NAM)+1:length(NAM)+4) = '.inp';
if STEP.Ana.Restart == 1; NAME.RestartName = NAM; NAME.RestartName(length(NAM)+1:length(NAM)+12) = '_restart.inp'; end

if GEO.cyl.Inc == 1
    hod = OP.Hc*0.5/OP.R; OP.Hod = hod; 
    if hod >= 2; OP.Class = 'Slender';
    elseif hod >= 1 && hod < 2; OP.Class = 'Intermediate';
    elseif hod >= 0.4 && hod < 1; OP.Class = 'Squat';
    else OP.Class = 'Retaining';
    end
else
    OP.Class = 'No silo';
end

OP.type = zeros(STEP.Tot,max(STEP.Sub));
OP.mag = zeros(STEP.Tot,max(STEP.Sub));
OP.hop_th0 = zeros(STEP.Tot,max(STEP.Sub));
OP.hop_eta = zeros(STEP.Tot,max(STEP.Sub));
OP.ears = zeros(STEP.Tot,max(STEP.Sub));
OP.factor = zeros(STEP.Tot,max(STEP.Sub));
OP.kc = zeros(STEP.Tot,max(STEP.Sub));
OP.chan0 = zeros(STEP.Tot,max(STEP.Sub));
for S = 1:STEP.Tot
    for U = 1:STEP.Sub(S)
        OP.type(S,U) = 0; OP.mag(S,U) = 0; OP.hop_th0(S,U) = 0; OP.hop_eta(S,U) = 0; OP.ears(S,U) = 0; OP.factor(S,U) = 0;
        for X = 1:size(STEP.Load.Type,1)
            if STEP.Load.Type(X,1:2) == [S U]; OP.type(S,U) = STEP.Load.Type(X,3); end
        end
        for X = 1:size(STEP.Load.A_mag,1)
            if STEP.Load.A_mag(X,1:2) == [S U]; OP.mag(S,U) = STEP.Load.A_mag(X,3); end
        end
        for X = 1:size(STEP.Load.Hopper_data,1)
            if STEP.Load.Hopper_data(X,1:2) == [S U]; OP.hop_th0(S,U) = STEP.Load.Hopper_data(X,3); OP.hop_eta(S,U) = STEP.Load.Hopper_data(X,4); end
        end
        for X = 1:size(STEP.Load.Ears,1)
            if STEP.Load.Ears(X,1:2) == [S U]; OP.ears(S,U) = STEP.Load.Ears(X,3); end
        end
        for X = 1:size(STEP.Load.Factor,1)
            if STEP.Load.Factor(X,1:2) == [S U]; OP.factor(S,U) = STEP.Load.Factor(X,3); end
        end
        for X = 1:size(STEP.Load.EN_Dis,1)
            if STEP.Load.EN_Dis(X,1:2) == [S U]; OP.EN_Dis(S,U) = STEP.Load.EN_Dis(X,3); end
        end
        for X = 1:size(STEP.Load.Channel_data,1)
            if STEP.Load.Channel_data(X,1:2) == [S U]; OP.chan0(S,U) = STEP.Load.Channel_data(X,3); OP.kc(S,U) = STEP.Load.Channel_data(X,4); end
        end
    end
end       

if ~isempty(find(STEP.Ana.Type == 3, 1))
    OP.Riks = 1;
else OP.Riks = 0;
end
if ~isempty(MATERIAL.sy); OP.plastic = 1; else OP.plastic = 0; end
OP.count_1 = 0; OP.count_2 = 0; OP.count_3 = 0;

OP.Heltoggle = 0;
if GEO.helix.Toggle == 1
    if GEO.helix.choice == 1
            N_p = OP.H/GEO.helix.pitch;
            Th_f = 2*pi*N_p;
            aa = OP.R; bb = GEO.helix.pitch/(2*pi);
            s = Th_f*sqrt(aa*aa+bb*bb);
            alpha = asin(OP.H/s);
            GEO.helix.incline = alpha*180/pi;
            GEO.helix.no_p = N_p;
    elseif GEO.helix.choice == 2
            alpha = GEO.helix.incline*pi/180;
            N_p = abs((0.5)*((-sin(alpha)^2+1)^0.5)*OP.H/(pi*OP.R*sin(alpha)));
            GEO.helix.pitch = OP.H/N_p;
            GEO.helix.no_p = N_p;
            Th_f = 2*pi*N_p;
            aa = OP.R; bb = GEO.helix.pitch/(2*pi);
            s = Th_f*sqrt(aa*aa+bb*bb);            
    end
    OP.Hpitch = GEO.helix.pitch; OP.N_p = GEO.helix.no_p;
    GEO.Ttot = 360;
    OP.Heltoggle = 1; OP.Helt = GEO.helix.thick;
    
    if GEO.helix.element == 'a' || GEO.helix.element == 'b' || GEO.helix.element == 'c'; % S4, S4R or S4R5
        GEO.helix.order = 1; OP.Helorder = GEO.helix.order;
    elseif GEO.helix.element == 'g' || GEO.helix.element == 'h' || GEO.helix.element == 'j'; % S8R, S8R5, S9R5
        GEO.helix.order = 2; OP.Helorder = GEO.helix.order;
    end
    OP.Heltype = GEO.helix.element; OP.T = GEO.Ttot; I = 0;
    if ~isempty(RSL.HEL.h1)
       for H = 1:length(RSL.HEL.h1); I = I + 1;
            RSL.Region(I) = I; RSL.Where(I) = 2; RSL.h1(I) = RSL.HEL.h1(H); RSL.h2(I) = RSL.HEL.h2(H);        
            RSL.t(I) = GEO.helix.thick; RSL.Type(I) = GEO.helix.element;            
            if GEO.helix.order == 1; RSL.Z(I) = RSL.HEL.Z(H); elseif GEO.helix.order == 2; RSL.Z(I) = RSL.HEL.Z(H)*2; end
       end
    end

    if GEO.helix.conformal == 0 
        for N = 1:ceil(GEO.helix.no_p+1)
        RSL.t1(N) = (N-1)*2*pi; RSL.t2(N) = N*2*pi; 
            if GEO.helix.order == 1
                RSL.T(N) = RSL.HEL.dp; 
            elseif GEO.helix.order == 2 
                RSL.T(N) = RSL.HEL.dp*2; 
            end        
        end
        RSL.Ttot = sum(RSL.T);
        RSL.h1 = RSL.h1*GEO.helix.pitch; RSL.h2 = RSL.h2*GEO.helix.pitch;
    elseif GEO.helix.conformal == 1
        arc = @(phi,P,R) phi*sqrt(R^2 + (P/(2*pi))^2);
        frac = (2*pi*OP.R*cos(GEO.helix.incline*pi/180))/arc(2*pi,GEO.helix.pitch,OP.R);  
        dp = RSL.HEL.dp;
        for N = 1:ceil(GEO.helix.no_p)           
            RSL.t1((2*N-1):(2*N)) = [0    frac]*arc(2*pi,GEO.helix.pitch,OP.R) + (N-1)*arc(2*pi,GEO.helix.pitch,OP.R);
            RSL.t2((2*N-1):(2*N)) = [frac 1   ]*arc(2*pi,GEO.helix.pitch,OP.R) + (N-1)*arc(2*pi,GEO.helix.pitch,OP.R);;
            RSL.T((2*N-1):(2*N)) = [round(frac*dp) round((1-frac)*dp)];
        end
        RSL.t1(end+1) = RSL.t2(end); RSL.t2(end+1) = RSL.t2(end) + RSL.t2(1); RSL.T(end+1) = RSL.T(1);
        if GEO.helix.order == 2; RSL.T = RSL.T*2; end        
        RSL.Ttot = sum(RSL.T);
        RSL.h1 = RSL.h1*2*pi*OP.R*sin(GEO.helix.incline*pi/180); RSL.h2 = RSL.h2*2*pi*OP.R*sin(GEO.helix.incline*pi/180);
    end
end