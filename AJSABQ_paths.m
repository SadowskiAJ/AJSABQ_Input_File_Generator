function [PATHS] = AJSABQ_paths(PATHS,NODES,GEO,OP,NAME,RSLTtot)
% Function to extract the paths for ABAQUS CAE data use 
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:24 (previously 22/11/09 - 00:46)

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


% Default values
dT = 45; % Degrees
dA = 0.25; % 25% of the height per component, except roof

disp(' '); disp('Extracting CAE paths.');
Hh = GEO.hop.Htot; Hc = GEO.cyl.Htot;
rangeZH = []; rangeZC = []; rangeT = []; H = 0; C = 0; T = 0;
rangeZHd = []; rangeZHa = []; rangeZCd = []; rangeZCa = [];
if GEO.hop.Inc == 1
    for a = 0:dA:1; H = H + 1;
        if a == 0; rangeZHd(H) = 1; else rangeZHd(H) = -1; end
        rangeZH(H) = a*Hh; rangeZHa(H) = a;
    end
    if a ~= 1; rangeZH(H+1) = 0.98*Hh; rangeZHd(H+1) = 21; rangeZHa(H+1) = 0.98; else rangeZHd(H) = 21; end
end
if GEO.cyl.Inc == 1
    for c = 0:dA:1; C = C + 1;
        if c == 0; rangeZCd(C) = 22; else rangeZCd(C) = -2; end
        rangeZC(C) = c*Hc + Hh; rangeZCa(C) = c;
    end
    if c ~= 1; rangeZC(C+1) = Hc + Hh; rangeZCd(C+1) = 3; rangeZCa(C+1) = 1; else rangeZCd(C) = 3; end
    if rangeZC(1) == 0 && GEO.hop.Inc == 1; rangeZC(1) = 0.02*Hc + Hh; rangeZCa(1) = 0.02; end
    if rangeZC(length(rangeZC)) == Hc && GEO.roof.Inc == 1;
        rangeZC(length(rangeZC)) = 0.99*Hc + Hh;
        rangeZCd(length(rangeZCd)) = 31; rangeZCa(length(rangeZCa)) = 0.99;
    end
end
rangeZ = [rangeZH rangeZC]; rangeZd = [rangeZHd rangeZCd]; rangeZa = [rangeZHa rangeZCa];

for t = 0:dT:GEO.Ttot; T = T + 1;
    rangeT(T) = t;
end
if t ~= GEO.Ttot; rangeT(T+1) = GEO.Ttot; end
if rangeT(1) == 0 && rangeT(length(rangeT)) == 360; rangeT(length(rangeT)) = []; end
PATHS.Circumferential = ''; 
fT = 0; TN = struct('I',[],'T',[],'Z',[]);
for T = 1:length(rangeT); fT = fT + 1;
    [dif in] = min(abs(NODES.T - rangeT(T)));
    I = find(NODES.T == NODES.T(in)); Ino = find(NODES.Z(I) > OP.Hc); I(Ino) = [];
    PATHS.Axial.Generatrix(T) = rangeT(T) + dif;
    PATHS.Axial.First(T) = NODES.index(I(1));
    PATHS.Axial.Last(T) = NODES.index(I(length(I)));
    PATHS.Axial.Step(T) = RSLTtot+1;
    if T == 1; PATHS.Axial.Description(T,:) = 'Left ';
    elseif T == length(rangeT); PATHS.Axial.Description(T,:) = 'Right';
    else PATHS.Axial.Description(T,:) = '     ';
    end
end

for Z = 1:length(rangeZ)
    [dif in] = min(abs(NODES.Z - rangeZ(Z)));
    I = find(NODES.Z == NODES.Z(in));
    PATHS.Circumferential.Generatrix(Z) = rangeZ(Z) + dif;
    PATHS.Circumferential.First(Z) = NODES.index(min(I));
    PATHS.Circumferential.Last(Z) = NODES.index(max(I));
    PATHS.Circumferential.Step(Z) = 1;
    PATHS.Circumferential.Fraction(Z,1) = rangeZa(Z); PATHS.Circumferential.Fraction(Z,2) = rangeZd(Z);
    switch rangeZd(Z)
        case 1
            PATHS.Circumferential.Description(Z,:) = ['Hopper outlet             '];
        case -1
            PATHS.Circumferential.Description(Z,:) = 'within hopper             ';
        case 21
            PATHS.Circumferential.Description(Z,:) = 'Transition (hopper side)  ';
        case 22
            PATHS.Circumferential.Description(Z,:) = 'Transition (cylinder side)';
        case -2
            PATHS.Circumferential.Description(Z,:) = 'within cylinder           ';
        case 3
            PATHS.Circumferential.Description(Z,:) = 'Top of silo               ';
        case 31
            PATHS.Circumferential.Description(Z,:) = 'Top of silo (below roof)  ';
    end
end

% If the ABAQUS CAE Python script has been requested
if PATHS.Python.Toggle == 1
    [PATHS] = AJSABQ_fwrite_python(PATHS,OP.Copy,NAME);
end
