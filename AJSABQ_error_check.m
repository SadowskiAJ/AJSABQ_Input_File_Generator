function [go] = AJSABQ_error_check(GEO,NAME,DRAW,MATERIAL,RSL,STEP,DATA,SOLID,SPRING,PATHS,PP)
% Function to check the consistency of the input data for the InputOP.type
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:24 (previously 28/02/12 - 17:08)

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

go = 1;
if GEO.cyl.Inc ~= 1 && GEO.hop.Inc ~= 1
    go = 0; disp('Either a hopper or a cylinder must be present.');
end
if GEO.hop.Inc == 1 && GEO.roof.Inc == 1 && GEO.cyl.Inc ~= 1
    go = 0; disp('A roof must go on a cylinder, not a hopper.');
end
if GEO.roof.Inc == 1
    if isempty(GEO.roof.Angle) || isempty(GEO.roof.Thick); go = 0; disp('No roof angle/thickness specified.'); end
    if isempty(RSL.CYL.RoofZ); go = 0; disp('No roof elements specified.'); end
end
if std([length(RSL.HOP.h1) length(RSL.HOP.h2) length(RSL.HOP.t) length(RSL.HOP.Type)]) ~= 0
    go = 0; disp('Hopper axial resolution struct vectors must have the same length.');
else
    if GEO.hop.Inc == 1
        for H = 1:length(RSL.HOP.Type)
            a = length(find(RSL.HOP.Type(H) == 'a')); b = length(find(RSL.HOP.Type(H) == 'b')); c = length(find(RSL.HOP.Type(H) == 'c'));
            d = length(find(RSL.HOP.Type(H) == 'd')); e = length(find(RSL.HOP.Type(H) == 'e')); f = length(find(RSL.HOP.Type(H) == 'f'));
            g = length(find(RSL.HOP.Type(H) == 'g')); h = length(find(RSL.HOP.Type(H) == 'h')); j = length(find(RSL.HOP.Type(H) == 'j'));
            if sum([a b c d e f g h j]) == 0; go = 0; disp('No such shell element (hopper).'); end
        end
    end
end
if std([length(RSL.CYL.h1) length(RSL.CYL.h2) length(RSL.CYL.t) length(RSL.CYL.Type)]) ~= 0
    go = 0; disp('Cylinder axial resolution struct vectors must have the same length.');
else
    if GEO.cyl.Inc == 1
        for H = 1:length(RSL.CYL.Type)
            a = length(find(RSL.CYL.Type(H) == 'a')); b = length(find(RSL.CYL.Type(H) == 'b')); c = length(find(RSL.CYL.Type(H) == 'c'));
            d = length(find(RSL.CYL.Type(H) == 'd')); e = length(find(RSL.CYL.Type(H) == 'e')); f = length(find(RSL.CYL.Type(H) == 'f'));
            g = length(find(RSL.CYL.Type(H) == 'g')); h = length(find(RSL.CYL.Type(H) == 'h')); j = length(find(RSL.CYL.Type(H) == 'j'));
            if sum([a b c d e f g h j]) == 0; go = 0; disp('No such shell element (cylinder).'); end
        end
    end
end
if std([length(RSL.CIRC.t1) length(RSL.CIRC.t2) length(RSL.CIRC.T)]) ~= 0
    go = 0; disp('Global circumferential resolution struct components must have the same length.');
end
if DATA.copy.Toggle == 1 && isempty(DATA.copy.Destination); go = 0; disp('Destination copy path for .inp file must be specified.'); end
if DATA.python.Toggle ~= 1 && DATA.python.Copy == 1; go = 0; disp('Python script must be requested if it is to be copied.'); end
if DATA.python.Toggle == 1 && DATA.python.Copy == 1 && isempty(DATA.copy.Destination); go = 0; disp('Destination copy path for .py script must be specified.'); end
if DATA.python.Toggle == 1 && DATA.paths.Toggle ~= 1; go = 0; disp('Path extraction must be requested in order to generate a Python script.'); end
if GEO.cyl.Inc ~= 1; len = length(RSL.CYL.h1);
    if ~isempty(find(find(STEP.Load.Type(:,3)) == 10) == 1); go = 0; disp('Edge load requested but cylinder not included.'); end
    if ~isempty(find(find(STEP.Load.Type(:,3)) == 15) == 1); go = 0; disp('Internal pressure requested but cylinder not included.'); end
    if ~isempty(find(find(STEP.Load.Type(:,3)) == 152) == 1); go = 0; disp('Internal friction requested but cylinder not included.'); end
    if isempty(RSL.CYL.h1) || isempty(RSL.CYL.h2) || isempty(RSL.CYL.Z) || isempty(RSL.CYL.t) || isempty(RSL.CYL.Type); go = 0; disp('Insufficient cylinder mesh data given.'); end
    if RSL.CYL.h2(length(RSL.CYL.h2)) ~= GEO.cyl.Htot; go = 0; disp('Cylinder resolution must be defined along its entire height.'); end
    if len > 1
        for H = 1:len-1
            if RSL.CYL.h1(H+1) ~= RSL.CYL.h2(H); go = 0; disp('Subsequent h1 entries must be equal to the previous h2 entries in the cylinder resolution struct.'); end
        end
    end
end
if GEO.hop.Inc ~= 1; len = length(RSL.HOP.h1);
    if ~isempty(find(find(STEP.Load.Type(:,3)) == 16) == 1); go = 0; disp('Internal pressure requested but hopper not included.'); end
    if ~isempty(find(find(STEP.Load.Type(:,3)) == 162) == 1); go = 0; disp('Internal friction requested but hopper not included.'); end
    if isempty(RSL.HOP.h1) || isempty(RSL.HOP.h2) || isempty(RSL.HOP.Z) || isempty(RSL.HOP.t) || isempty(RSL.HOP.Type); go = 0; disp('Insufficient hopper mesh data given.'); end
    if RSL.HOP.h2(length(RSL.HOP.h2)) ~= GEO.hop.Htot; go = 0; disp('Hopper resolution must be defined along its entire height.'); end
    if len > 1
        for H = 1:len-1
            if RSL.HOP.h1(H+1) ~= RSL.HOP.h2(H); go = 0; disp('Subsequent h1 entries must be equal to the previous h2 entries in the hopper resolution struct.'); end
        end
    end
end
len = length(RSL.CIRC.t1); if length(RSL.CIRC.t2) ~= len || length(RSL.CIRC.T) ~= len; go = 0; disp('Circumferential mesh resolution vectors must be of the same length.'); end
if len > 1
    for H = 1:len-1
        if RSL.CIRC.t1(H+1) ~= RSL.CIRC.t2(H); go = 0; disp('Subsequent t1 entries must be equal to the previous t2 entries in the circumferential resolution struct.'); end
    end
end
if GEO.feature.C_Weld.Toggle == 1; len = length(GEO.feature.C_Weld.Type);
    if length(GEO.feature.C_Weld.Z) ~= len || length(GEO.feature.C_Weld.Ampt) ~= len; go = 0; disp('Incorrect circumferential weld data.'); end
    if length(find(GEO.feature.C_Weld.Type == 'A')) + length(find(GEO.feature.C_Weld.Type == 'B')) ~= len; go = 0; disp('Incorrect circumferential weld definition.'); end
    if GEO.feature.C_Weld.Z(len) > GEO.cyl.Htot; go = 0; disp('Circumferential welds extend beyond the cylinder height.'); end
end
if GEO.feature.A_Weld.Toggle == 1; len = length(GEO.feature.A_Weld.Z); RSLt = [RSL.HOP.t RSL.CYL.t];
    if length(GEO.feature.A_Weld.Ampt) ~= length(GEO.feature.A_Weld.t); go = 0; disp('Incorrect vertical weld data.'); end
    if GEO.feature.A_Weld.Type == 'A' || GEO.feature.A_Weld.Type == 'B'; else; go = 0; disp('Incorrect vertical weld definition.'); end
    if GEO.feature.A_Weld.Z(len) > GEO.cyl.Htot; go = 0; disp('Vertical welds extend beyond the cylinder height.'); end
    for A = 1:length(GEO.feature.A_Weld.t); goo = 0;
        for R = 1:length(RSLt)
            if RSLt(R) == GEO.feature.A_Weld.t(A); goo = 1; end
        end
        if goo == 0; go = 0; disp('Not all possible shell wall thicknesses have been accounted for in the vertical weld definition.'); break; end
    end
end
if ~isempty(GEO.feature.Super_ell) && length(GEO.feature.Super_ell) ~= 7; go = 0; disp('Insufficient super-elliptical flattening definition.'); end
if GEO.feature.Super_ell(1) == 1 && length(GEO.feature.Super_ell) == 7
    if GEO.Rtot*(1 - cos(GEO.feature.Super_ell(3)*pi/180)) < GEO.feature.Super_ell(2); go = 0; disp('Super-elliptical flattening depression/spread misalignment.'); end
    if GEO.feature.Circ_flat(1) == 1; go = 0; disp('Circular and super-elliptical flattening are mutually exclusive.'); end
    if GEO.feature.Super_ell(6) >= GEO.cyl.Htot || GEO.feature.Super_ell(6) <= 0; go = 0; disp('Position of max. axial imperfection out of bounds (super-elliptical flattening.'); end
end
if GEO.feature.Circ_flat(1) == 1
    Tc = GEO.feature.Circ_flat(3)*pi/180; d = GEO.feature.Circ_flat(2); R = GEO.Rtot;
    D = -(1/2)*d*(2*R-d)/(-R+d+R*cos(Tc));
    Rp = (1/2)*(-2*R^2+2*R*d-d^2+2*cos(Tc)*R^2-2*cos(Tc)*R*d)/(-R+d+R*cos(Tc));
    Tcp = (180/pi)*atan2(2*R*sin(Tc)*(-R+d+R*cos(Tc))/(-2*R^2+2*R*d-d^2+2*cos(Tc)*R^2-2*cos(Tc)*R*d), (-2*R*d+2*R^2*cos(Tc)^2+d^2+2*cos(Tc)*R*d-2*cos(Tc)*R^2)/(-2*R^2+2*R*d-d^2+2*cos(Tc)*R^2-2*cos(Tc)*R*d));
    if d >= R*(1 - cos(Tc)); go = 0; disp('Circular flattening to deep.'); end
    if GEO.feature.Super_ell(1) == 1; go = 0; disp('Circular and super-elliptical flattening are mutually exclusive.'); end
    if GEO.feature.Circ_flat(5) >= GEO.cyl.Htot || GEO.feature.Circ_flat(5) <= 0; go = 0; disp('Position of max. axial imperfection out of bounds (circular flattening).'); end
end
if ~isempty(GEO.feature.Axi_sin) && length(GEO.feature.Axi_sin) ~= 3; go = 0; disp('Insufficient inward axial sine wave definition.'); end
if GEO.feature.Axi_sin(1) == 1 && length(GEO.feature.Axi_sin) == 3
    if GEO.feature.Axi_sin(2) >= GEO.Rtot; go = 0; disp('Axial sine wave too deep.'); end
    if GEO.feature.Axi_sin(3) <= 0; go = 0; disp('Half-wave number of axial sine wave must be >= 1.'); end
end
if ~isempty(GEO.feature.Circ_cos) && length(GEO.feature.Circ_cos) ~= 4; go = 0; disp('Insufficient circumferential cosine wave definition.'); end
if GEO.feature.Circ_cos(1) == 1 && length(GEO.feature.Circ_cos) == 4
    if GEO.feature.Circ_cos(2) >= GEO.Rtot; go = 0; disp('Circumferential cosine wave too deep.'); end
    if GEO.feature.Circ_cos(3) <= 0; go = 0; disp('Wave number of circumferential cosine wave must be >= 1.'); end
end
if isempty(DRAW.line_2D) ~= 1
    for I = 1:size(DRAW.line_2D,1); if DRAW.line_2D(I,1) > length(STEP.Ana.Type); go = 0; disp('2D line plot requested for a step which has not been defined.'); end; end
    if size(DRAW.line_2D,2) ~= 2; go = 0; disp('Incorrect 2D line plot request matrix size.'); end
end
if isempty(DRAW.surf_3D.step) ~= 1
    for I = 1:length(DRAW.surf_3D.step)
        if DRAW.surf_3D.step(I) > length(STEP.Ana.Type); go = 0; disp('3D surface plot requested for a step which has not been defined.'); end
    end
    if length(DRAW.surf_3D.which) ~= length(DRAW.surf_3D.step); go = 0; disp('3D surface .which vector must be the same length as the .step vector');
    else
        for I = 1:length(DRAW.surf_3D.step)
            if DRAW.surf_3D.which(I) ~= 1 && DRAW.surf_3D.which(I) ~= 2 && DRAW.surf_3D.which(I) ~= 3; go = 0; disp('3D surface incorrect request'); end
            pos = find(DRAW.surf_3D.step == DRAW.surf_3D.step(I));
            if length(pos) > 1; go = 0; disp('3D surface step requests cannot be duplicated.'); end
        end
    end
end
if isempty(DRAW.planar_2D.step) ~= 1
    for I = 1:length(DRAW.planar_2D.step)
        if DRAW.planar_2D.step(I) > length(STEP.Ana.Type); go = 0; disp('2D planar plot requested for a step which has not been defined.'); end
    end
    for I = 1:length(DRAW.planar_2D.step)
        pos = find(DRAW.planar_2D.step == DRAW.planar_2D.step(I));
        if length(pos) > 1; go = 0; disp('2D planar step requests cannot be duplicated.'); end
    end
end
for S = 1:length(STEP.Ana.Type)
    if isempty(find(STEP.Load.Type(:,1) == S, 1)); go = 0; disp('There must be at least one load in every step.'); end
end
if isempty(STEP.Ana.GN); go = 0; disp('Geometric (non)linearity not specified.'); end
if isempty(NAME.FileName); go = 0;
    disp('Please enter an output file name.');
end
if length(MATERIAL.sy) ~= length(MATERIAL.ep); go = 0;
    disp('Verify tabular material data.');
end
if ~isempty(MATERIAL.sy) && ~isempty(MATERIAL.ep) && MATERIAL.ep(1) ~= 0; go = 0;
    disp('Initial plastic strains must be zero.');
end
if STEP.Ana.Restart == 1 && (length(STEP.Ana.Riks_stepR) ~= length(STEP.Ana.Riks_incR)); go = 0;
    disp('Verify restart Riks data.');
end
if STEP.Ana.Restart == 1
    if STEP.Ana.Type(end) ~= 3 || STEP.Ana.GN(end) ~= 1; go = 0;
        disp('Restart file can only be requested if the final previous step is a nonlinear Riks step.');
    end
end
if ~isempty(find(STEP.Ana.Type == 2, 1)); loc = find(STEP.Ana.Type == 2);
    if length(STEP.Ana.LBA_no) < 2; go = 0;
        disp('Insufficient perturbation data');
    else
        for L = 1:length(loc)
            if isempty(find(STEP.Ana.LBA_no(:,1) == loc(L), 1)); go = 0;
                disp('Insufficient perturbation data.');
            end
        end
    end
end
if ~isempty(find(STEP.Ana.Type == 3, 1)); loc = find(STEP.Ana.Type == 3);
    for L = 1:length(loc)
        if isempty(find(STEP.Ana.Riks_step(:,1) == loc(L), 1)) || isempty(find(STEP.Ana.Riks_inc(:,1) == loc(L), 1)); go = 0;
            disp('Insufficient Riks data.');
        end
        if isempty(find(STEP.Ana.Riks_step(:,1) == loc(L), 1)) || isempty(find(STEP.Ana.Riks_max(:,1) == loc(L), 1)); go = 0;
            disp('Insufficient Riks data.');
        end
    end
end
if length(find(STEP.Ana.Type == 1)) + length(find(STEP.Ana.Type == 2)) == length(STEP.Ana.Type)
    if STEP.Ana.Restart == 1; go = 0;
        disp('Restart file requested for non-Riks step.');
    end
end

if size(STEP.Load.Type,1) == 1
    if STEP.Load.Type(1) ~= 1 || STEP.Load.Type(2) ~= 1; go = 0; disp('Incorrect load assignment.'); end
end
if size(STEP.Load.Channel_data,2) ~= 4 && size(STEP.Load.Channel_data,2) ~= 0; go = 0; disp('Incorrect channel data size.'); end
if size(STEP.Load.Type,1) > 0
    cur = 1; lds = []; no = 0;
    for L = 1:size(STEP.Load.Type,1)
        if STEP.Load.Type(L,1) == cur; lds(length(lds)+1) = STEP.Load.Type(L,3); end
        if STEP.Load.Type(L,1) ~= cur || L == size(STEP.Load.Type,1); ldsO = lds; lds = []; lds(length(lds)+1) = STEP.Load.Type(L,3);
            cur = STEP.Load.Type(L,1);
            for L1 = 1:length(ldsO); ldsOm = ldsO; ldsOm(L1) = [];
                for L2 = 1:length(ldsOm)
                    if ldsO(L1) == ldsOm(L2); go = 0; disp('Load repetition.'); end
                end
            end
            if length(find([length(find(ldsO == 19)) length(find(ldsO == 20)) length(find(ldsO == 201)) length(find(ldsO == 25)) length(find(ldsO == 30)) length(find(ldsO == 35)) length(find(ldsO == 50)) length(find(ldsO == 51)) length(find(ldsO == 100))] == 0)) < 8; go = 0;
                disp('Mutually exclusive silo load assignments.');
            end
            if length(find([length(find(ldsO == 16)) length(find(ldsO == 162)) length(find(ldsO == 40)) length(find(ldsO == 45)) length(find(ldsO == 46))] == 0)) < 4; go = 0;
                disp('Mutually exclusive hopper load assignments.');
            end
        end

        if STEP.Load.Type(L,3) == 10 || STEP.Load.Type(L,3) == 15 || STEP.Load.Type(L,3) == 152 || STEP.Load.Type(L,3) == 16 || STEP.Load.Type(L,3) == 162 ||...
                STEP.Load.Type(L,3) == 19 || STEP.Load.Type(L,3) == 20 || STEP.Load.Type(L,3) == 201 || STEP.Load.Type(L,3) == 25 || STEP.Load.Type(L,3) == 30 || STEP.Load.Type(L,3) == 35
            for X = 1:size(STEP.Load.Channel_data,1)
                if STEP.Load.Channel_data(X,1:2) == STEP.Load.Type(L,1:2); no = no + 1; end
            end
        end
        %if L == size(STEP.Load.Type,1) && no ~= size(STEP.Load.Channel_data,1); go = 0; disp('Incorrect channel data.'); end

        no = 0;
        if STEP.Load.Type(L,3) == 40 || STEP.Load.Type(L,3) == 45 || STEP.Load.Type(L,3) == 46
            for X = 1:size(STEP.Load.Hopper_data,1)
                if STEP.Load.Hopper_data(X,1:2) == STEP.Load.Type(L,1:2); no = no + 1; end
            end
        end
        if L == size(STEP.Load.Type,1) && no ~= size(STEP.Load.Hopper_data,1) && GEO.hop.Inc == 1; go = 0; disp('Incorrect hopper data.'); end

        if STEP.Load.Type(L,3) == 10 || STEP.Load.Type(L,3) == 15 || STEP.Load.Type(L,3) == 152 || STEP.Load.Type(L,3) == 16 || STEP.Load.Type(L,3) == 162; p = 0;
            for X = 1:size(STEP.Load.A_mag,1)
                if STEP.Load.A_mag(X,1:2) == STEP.Load.Type(L,1:2); p = 1; end
            end
            if p == 0; go = 0; disp('Load magnitudes not defined.'); end
        end

        if STEP.Load.Type(L,3) == 19 || STEP.Load.Type(L,3) == 25 || STEP.Load.Type(L,3) == 35; p = 0;
            for X = 1:size(STEP.Load.Channel_data,1)
                if STEP.Load.Channel_data(X,1:2) == STEP.Load.Type(L,1:2); p = 1; end
            end
            if p == 0; go = 0; disp('Silo channel data not defined.'); end
        end

        if STEP.Load.Type(L,3) == 45 || STEP.Load.Type(L,3) == 46; p = 0;
            for X = 1:size(STEP.Load.Hopper_data,1)
                if STEP.Load.Hopper_data(X,1:2) == STEP.Load.Type(L,1:2); p = 1; end
            end
            if p == 0; go = 0; disp('Hopper data not defined.'); end
        end

        if L == 1
            if STEP.Load.Type(L,1) ~= 1 || STEP.Load.Type(L,2) ~= 1; go = 0; disp('Incorrect load assignment.'); end
        else
            if (STEP.Load.Type(L,1) == STEP.Load.Type(L-1,1)) && (STEP.Load.Type(L,2) ~= STEP.Load.Type(L-1,2) + 1)
                go = 0; disp('Incorrect load assignment.');
            elseif (STEP.Load.Type(L,1) ~= STEP.Load.Type(L-1,1)) && (STEP.Load.Type(L,1) ~= STEP.Load.Type(L-1,1) + 1)
                go = 0; disp('Incorrect load assignment.');
            elseif (STEP.Load.Type(L,1) ~= STEP.Load.Type(L-1,1)) && STEP.Load.Type(L,2) ~= 1
                go = 0; disp('Incorrect load assignment.');
            end
        end
    end
end

if size(find(STEP.Load.A_mag(:,3) == 1),1) == 1 || size(find(STEP.Load.A_mag(:,3) == 2),1) == 1 || size(find(STEP.Load.A_mag(:,3) == 3),1) == 1
    if isempty(STEP.BC.Top_RigidRing); go = 0; disp('Point loads through rigid point not allowed if rigit point not requested (top).'); end
end
if size(find(STEP.Load.A_mag(:,3) == 4),1) == 1 || size(find(STEP.Load.A_mag(:,3) == 5),1) == 1 || size(find(STEP.Load.A_mag(:,3) == 6),1) == 1
    if isempty(STEP.BC.Top_RigidRing); go = 0; disp('Point moments through rigid point not allowed if rigit point not requested (top).'); end
end
if size(find(STEP.Load.A_mag(:,3) == 7),1) == 1 || size(find(STEP.Load.A_mag(:,3) == 8),1) == 1 || size(find(STEP.Load.A_mag(:,3) == 9),1) == 1
    if isempty(STEP.BC.Top_RigidRing); go = 0; disp('Point rotations through rigid point not allowed if rigit point not requested (top).'); end
end

if size(find(STEP.Load.A_mag(:,3) == .1),1) == 1 || size(find(STEP.Load.A_mag(:,3) == .2),1) == 1 || size(find(STEP.Load.A_mag(:,3) == .3),1) == 1
    if isempty(STEP.BC.Trans_RigidRing); go = 0; disp('Point loads through rigid point not allowed if rigit point not requested (bottom).'); end
end
if size(find(STEP.Load.A_mag(:,3) == .4),1) == 1 || size(find(STEP.Load.A_mag(:,3) == .5),1) == 1 || size(find(STEP.Load.A_mag(:,3) == .6),1) == 1
    if isempty(STEP.BC.Trans_RigidRing); go = 0; disp('Point moments through rigid point not allowed if rigit point not requested (bottom).'); end
end
if size(find(STEP.Load.A_mag(:,3) == .7),1) == 1 || size(find(STEP.Load.A_mag(:,3) == .8),1) == 1 || size(find(STEP.Load.A_mag(:,3) == .9),1) == 1
    if isempty(STEP.BC.Trans_RigidRing); go = 0; disp('Point rotations through rigid point not allowed if rigit point not requested (bottom).'); end
end

if STEP.Load.Type(L,3) == 19 || STEP.Load.Type(L,3) == 20 || STEP.Load.Type(L,3) == 201 || STEP.Load.Type(L,3) == 25 || STEP.Load.Type(L,3) == 30 || STEP.Load.Type(L,3) == 35
    if sum(find(SOLID.Wall == 'D1')) + sum(find(SOLID.Wall == 'D2')) + sum(find(SOLID.Wall == 'D3')) ~= 5; go = 0; disp('Incorrect wall smoothness assignment (D1, D2 or D3 only).'); end
end
if ~isempty(STEP.Load.Ears) && size(STEP.Load.Ears,2) ~= 3
    go = 0; disp('Incorrect input (Ears)');
end
if ~isempty(STEP.Load.Factor) && size(STEP.Load.Factor,2) ~= 3
    go = 0; disp('Incorrect input (Factor)');
end
if ~isempty(STEP.Load.EN_Dis) && size(STEP.Load.EN_Dis,2) ~= 3
    go = 0; disp('Incorrect input (EN Discharge factors)');
end
if ~isempty(find(STEP.Load.Type(:,3) == 100, 1)) && isempty(STEP.Load.QVV)
    go = 0; disp('External .qvv file not defined');
end

if GEO.helix.Toggle ~= 1
    if size(STEP.BC.T0,1) ~= length(STEP.Ana.Type)
        if max(RSL.CIRC.t2) ~= 360; go = 0; disp('Incorrect boundary conditions at T = 0.'); end
    end
    if size(STEP.BC.Ttot,1) ~= length(STEP.Ana.Type)
        if max(RSL.CIRC.t2) ~= 360; go = 0; disp('Incorrect boundary conditions at T = Ttot.'); end
    end    
end
if GEO.hop.Inc ~= 1 && isempty(STEP.BC.Pit) ~= 1
    go = 0; disp('Hopper not defined but BC for hopper requested.');
end
if MATERIAL.Type == 2
    if length(MATERIAL.E) ~= 2 || length(MATERIAL.G) < 1; go = 0; disp('Incorrect material moduli definition.'); end
end
if length(MATERIAL.sy) ~= length(MATERIAL.ep); go = 0; disp('Incorrect material plasticity data definition.'); end
if length(find(STEP.Ana.Type == 1)) + length(find(STEP.Ana.Type == 2)) + length(find(STEP.Ana.Type == 3)) == 0; go = 0; disp('Incorrect step type definition.'); end
if SPRING.toggle == 1 && SPRING.Discharge == 1
    if std(STEP.Load.Type(:,3)) > 0; go = 0; disp('Only one discharge load case must be present with springs.'); end
    if std(STEP.Load.Channel_data(:,3)) > 0; go = 0; disp('Only one discharge load case must be present with spring.'); end
    if length(find(STEP.Load.Type == 10)) + length(find(STEP.Load.Type == 15)) + length(find(STEP.Load.Type == 152)) > 0; go = 0; disp('Discharge springs not defined for this case (cylinder).'); end
    if length(find(STEP.Load.Type == 16)) + length(find(STEP.Load.Type == 162)) + length(find(STEP.Load.Type == 40)) +...
            length(find(STEP.Load.Type == 45)) + length(find(STEP.Load.Type == 46)) > 0; go = 0; disp('Discharge springs not defined for this case (hopper).'); end
end
if PP.CPUs > str2num(getenv('NUMBER_OF_PROCESSORS'))
    go = 0; disp('Too many processors requested.');
end

if ~isempty(STEP.BC.Top_RigidRing) && GEO.roof.Inc == 1; go = 0; disp('Cannot request roof and top rigid ring together.'); end
if ~isempty(STEP.BC.Top_RigidRing) && ~isempty(STEP.BC.Top); go = 0; disp('Only one of top shell BC or top rigid ring BC must be defined.'); end
if ~isempty(STEP.BC.Trans_RigidRing) && GEO.hop.Inc == 1; go = 0; disp('Cannot request hopper and bottom rigid ring together.'); end
if ~isempty(STEP.BC.Trans_RigidRing) && ~isempty(STEP.BC.Trans); go = 0; disp('Only one of bottom shell BC or bottom rigid ring BC must be defined.'); end

% for L = 1:size(DATA.Load_Type,1);
%     if DATA.Load_Type(L,3) == 17; % i.e. if it's the inline function
%         if L > 1 && DATA.Load_Type(L-1,1) == DATA.Load_Type(L,1); go = 0; disp('Custom inline function must be defined as the only load case.'); end
%         if L < size(DATA.Load_Type,1) && DATA.Load_Type(L+1,1) == DATA.Load_Type(L,1); go = 0; disp('Custom inline function must be defined as the only load case.'); end
%     end
% end

if ~isempty(find(DATA.Load_Type(:,3) == 201, 1)) && isempty(DATA.Load.Inline.F); go = 0; disp('Please define a custom inline load function.'); end
el = 0;
if GEO.helix.Toggle == 1
    len = length(RSL.HEL.h1);
    if len > 1
        if RSL.HEL.h1(1) ~= 0; go = 0; disp('Helical axial resolution struct is normalised between 0 and 1. Starts at 0.'); end
        if RSL.HEL.h2(len) ~= 1; go = 0; disp('Helical axial resolution struct is normalised between 0 and 1. Ends at 1.'); end
    end
    if mod(RSL.HEL.Z,1) > 0; go = 0; disp('Number of requested helical (''vertical'') elements must be an integer.'); end
    if mod(RSL.HEL.dp,1) > 0; go = 0; disp('Number of requested helical (''arc length'') elements must be an integer.'); end
    if GEO.roof.Inc == 1; go = 0; disp('Cannot request a roof with helical meshing.'); end
    if GEO.hop.Inc == 1; go = 0; disp('Cannot request a hopper with helical meshing.'); end
    if SPRING.toggle == 1; go = 0; disp('Cannot request spring elements with helical meshing.'); end
    
    if sum(size(STEP.BC.T0)) > 0 || sum(size(STEP.BC.Ttot)) > 0
        go = 0; disp('Helical meshing requires the full shell. Don''t specify circumferential BCs.'); 
    end
    
    if GEO.feature.C_Weld.Toggle == 1
        go = 0; disp('Do not request circumferential welds with helical meshing.'); 
    end
    if GEO.feature.A_Weld.Toggle == 1 
        go = 0; disp('Do not request axial welds with helical meshing.'); 
    end
    if GEO.feature.H_Weld.Toggle == 1
        if length(find(GEO.feature.H_Weld.Type == 'A')) + length(find(GEO.feature.H_Weld.Type == 'B')) == 0
            go = 0; disp('Only type A or B spiral welds accepted.');
        end            
    end
    if GEO.helix.element == 'a'; el = 1; end
    if GEO.helix.element == 'b'; el = 1; end
    if GEO.helix.element == 'c'; el = 1; end
    if GEO.helix.element == 'g'; el = 1; end
    if GEO.helix.element == 'h'; el = 1; end
    if GEO.helix.element == 'j'; el = 1; end
    if el == 0; go = 0; disp('No appropriate helical element has been specified.'); end
    if PATHS.Toggle == 1; go = 0; disp('Cannot request nodal paths with helical meshing.'); end
    if sum(RSL.HEL.Z)/PP.CPUs < 1; go = 0; disp('Helical mesh is too small for the requested no. of CPUs.'); end
end

if go == 1; disp('      Input verified.'); end