function [GEO,NAME,DRAW,MATERIAL,RSL,STEP,SOLID,PATHS,SPRING] = AJSABQ_structs()
% Function to load the global struct definitions 
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:27 (previously 11/02/12 - 15:03)

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

GEO = struct('cyl',[],'hop',[],'roof',[],'helix',[],'feature',[],'Rtot',[],'Ttot',[]); % Global Geometry Struct
GEO.cyl = struct('Htot',[],'Inc',[]); % Vertical height and total radius; Include (1 - yes)
GEO.hop = struct('Htot',[],'Beta',[],'Inc',[]); % Vertical height and apex half-angle; Include (1 - yes)
GEO.roof = struct('Angle',[],'Thick',[],'Inc',[]); % Angle to the horizontal and thickness of roof; Include (1 - yes)
GEO.helix = struct('Toggle',[],'choice',[],'pitch',[],'incline',[],'no_p',[],'order',[],'conformal',[],'element',[]); % Request info on helical meshing
GEO.feature = struct('LBA',[],'H_Weld',[],'C_Weld',[],'A_Weld',[],'Super_ell',[],'Circ_flat',[],'Axi_sin',[],'Circ_cos',[],'Plate_level',[],'Cyl_perturbation',[]); % Secondary geometrical features
GEO.feature.H_Weld = struct('Toggle',[],'Type',[],'Ampt',[]); % Data of the Spiral Weld
GEO.feature.C_Weld = struct('Toggle',[],'Type',[],'Z',[],'Ampt',[]); % Data of the Circumferential Weld
GEO.feature.A_Weld = struct('Toggle',[],'Type',[],'Z',[],'Ampt',[],'t',[],'Tfreq',[],'Tphase',[]); % Data of the Axial Weld
NAME = struct('Heading',[],'FileName',[],'RestartName',[],'Description',[]); % File Name Struct
DRAW = struct('line_2D',[],'surf_3D',[],'mesh_3D',[]); % Drawing Request Struct
DRAW.planar_2D = struct('step',[]); % 2D planar detailed Request Struct
DRAW.surf_3D = struct('step',[],'which',[]); % 3D surf detailed Request Struct
DRAW.mesh_3D = struct('toggle',[],'nodes',[],'elements',[]); % 3D mesh detailed Request Struct
MATERIAL = struct('Name',[],'Type',[],'E',[],'v',[],'G',[],'sy',[],'ep',[]); % Material Definition Struct
RSL = struct('HOP',[],'CYL',[],'HEL',[],'CIRC',[],'Region',[],'Where',[],'h1',[],'h2',[],'Z',[],'t',[],'Type',[],'Els',[],'Etot',[],'Hoptot',[],'Cyltot',[],'Rooftot',[],'t1',[],'t2',[],'T',[],'Ttot',[],'RoofZ',[]); % Mesh Resolution Struct
RSL.HOP = struct('h1',[],'h2',[],'Z',[],'t',[],'Type',[],'RoofZ',[]); % Axial Mesh Resolution Struct for the Hopper (and for the roof)
RSL.CYL = struct('h1',[],'h2',[],'Z',[],'t',[],'Type',[]); % Axial Mesh Resolution Struct for the Cylinder
RSL.HEL = struct('h1',[],'h2',[],'Z',[],'dp',[]); % Mesh Resolution Struct for a Cylinder with helical meshing
RSL.CIRC = struct('t1',[],'t2',[],'T',[]); % Circumferential Mesh Resolution Struct for the entire structure
STEP = struct('Tot',[],'Sub',[],'Load',[],'BC',[],'Output',[]); % Step Definition Struct
STEP.Ana = struct('Type',[],'GN',[],'LBA_no',[],'Riks_step',[],'Riks_inc',[],'Restart',[],'Riks_stepR',[],'Riks_incR',[]); % Step no. ID
STEP.Load = struct('Type',[],'A_mag',[],'Ears',[],'Hopper_data',[],'Channel_data',[],'Factor',[],'QVV',[]); % Step no. ID & inside the step ID
STEP.BC = struct('Top',[],'Trans',[],'T0',[],'Ttot',[],'Pit',[]);
SOLID = struct('Name',[],'Wall',[],'Weight',[],'K',[],'Repose',[],'Equiv',[],'mewU',[],'mewL',[],'Ch',[],'Cw',[],'Frang',[]); % Solid Definition Struct (if needed)
PATHS = struct('Toggle',[],'Circumferential',[],'Axial',[],'Python',[]); % ABAQUS data extraction paths Struct
PATHS.Circumferential = struct('Generatrix',[],'Description','','First',[],'Last',[],'Step',[],'Fraction',[]); % For the circumferential paths (every 45 degrees)
PATHS.Axial = struct('Generatrix',[],'Description','','First',[],'Last',[],'Step',[]); % For the axial paths (every 0.25 of height component)
PATHS.Python = struct('Toggle',[],'Copy',[],'Result',[]); % To create a Python script for CAE and to copy it to the relevant directory
SPRING = struct('toggle',[],'index',[],'nodes',[],'region',[],'set',[],'constant',[],'discharge',[]); % Struct listing spring elements