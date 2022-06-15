function [DATA] = AJSABQ_input(); DATA.go = 1;
% Function to hold the input data 
% AJSABQ Input File Generator Version 2.0

% SPECIFICALLY this file sets up a helically meshed cylinder under axial
% compression

% To be run from the Matlab workspace command:
% clear; MatMPI_Delete_all; eval(MPI_Run('AJSABQ_drivrutin',N,{}))
% where N is the no. of threads to be created (usually equal to the cores on the CPU)
% NOTE this requires that this file is called AJSABQ_input.m

% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:45 (previously 01/07/12 - 21:04)

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

F_load = 0; % Load a dataset from the .mat file? 1 for yes, anything else will skip to line 24
F_save = 0; % Save the data to a .mat file? if so, specify name (only for F_load ~= 1)
F_run = 1;

if F_load == 1
    dir *.mat; F_file = input('Choose .mat file, type file name only with no extension. Type 0 to exit.');
    if F_file ~= 0 
        F_file(length(F_file)+1:length(F_file)+4) = '.mat'; load F_file DATA;
    elseif F_file == 0
        DATA.go = 0; disp('Load .mat file not chosen.');
    end
end
% Materials
E = 2e5; v = 0.3;

% Helix
angle = 45;

% Geometry
t = 1; R = 1000*t; D = 2*R; turns = 3;
%ZZ = 1000; L = ceil(sqrt(ZZ*R*t/(sqrt(1-v*v))));
L = round(turns*pi*tan(angle*pi/180)*D);
ampt = 3.0; 

% Data
scr = E*t/(R*sqrt(3*(1-v*v)));
Pcr = 2*pi*R*t*scr;
Ncr = t*scr;

% Analysis
stp = 0.02; 

if F_load ~= 1
    DATA.cInc = 1; % Include cylinder (1 - yes)
    DATA.hInc = 0; % Include hopper (1 - yes)
    DATA.rInc = 0; % Include roof (1 - yes)
    %DATA.Name = 'test'; % .inp file name
    DATA.Name = 'AJS_GNIA_Conf_45_d0t_3p00'; % .inp file name
    DATA.Heading = 'Helical meshing - NonConformal - R&T Type A Spiral Weld - S9R5 - 3 turns 15 degrees - R/t = 1000'; % ABAQUS description
    DATA.Description = 'Axial Compression GNIA'; % private description
    DATA.Rtot = R; % Total max. radius
    DATA.cHtot = L; % Height of cylinder (if required)
    DATA.hHtot = 0; % Height of hopper (if required)
    DATA.rAngle = 1e-6; % Slope to horizontal of roof (if required - degrees)
    DATA.rThick = 6; % Thickness of roof

    DATA.line_2D = []; % Draw a 2D line plot? [Step  (1 - yes)]
    DATA.planar_2D.step = []; % Draw a 2D planar plot of the flow channel geometry (only with appropriate load cases)? [Give vector of step nos., not exceeding the amount you have defined]
    DATA.surf_3D.step = [1]; % Draw a 3D surface plot? [Give vector of step nos., not exceeding the amount you have defined]
    DATA.surf_3D.which = [1]; % As long as the above step vectors are defined; 1 - normal pressures only; 2 - frictional tractions only; 3 - both separately
    DATA.mesh_3D.toggle = 0; % Draw a 3D line plot of the FEA mesh? (1 - yes)
    DATA.mesh_3D.nodes = 0; % Number nodes? (1 - yes); WARNING! MAY CRASH FOR LARGE MESHES
    DATA.mesh_3D.elements = 0; % Number elements? (1 - yes); WARNING! MAY CRASH FOR LARGE MESHES
    DATA.copy.Toggle = 1; % Copy the final .inp file to a destination?
    DATA.copy.Destination = 'S:\jobs'; % Full path name of destination, as a string, usually: 'F:\MATLAB7\work'
    DATA.paths.Toggle = 0; % Requests circumferential & axial paths for ABAQUS CAE data extraction? (every 45 degrees & 0.25H)
    DATA.python.Toggle = 0; % Requests that these be implemented in a Python script (paths MUST be requested for this to be generated)
    DATA.python.Copy = 0; % Copy the final .py file to the same destination as the .inp file?

    DATA.helical.Toggle = 1; % Request helical meshing
    DATA.helical.choice = 2; % Helix specification (1 or 2)
    DATA.helical.pitch = DATA.cHtot/3; % CHOICE 1 = Specify the pitch directly (height drop every 2Pi radians)
    DATA.helical.incline = angle; % CHOICE 2 = Specify the mesh inclination angle (degrees) to the horizontal
    DATA.helical.thick = t; % Uniform shell thickness
    DATA.helical.conformal = 1; % Generate conformal mesh?   
    DATA.helical.element = 'j'; % Dominant element type (see list below)
    %% NOTE1: choosing S4 ('a'), S4R ('b') or S4R5 ('c') will, of course, generate a 1st order mesh, with S3 elements as a filler
    %% NOTE2: choosing S8R ('g'), S8R5 ('h') or S9R5 ('j') will, of course, generate a 2nd order mesh, with STRI65 elements as a filler
    %% NOTE3: choosing S3 ('d'), STRI3 ('e') or STRI65 ('f') explicitely will do nothing
    
    % Helical meshing struct - generates full 360 degrees
    DATA.HELh1 = [0  ]; % For a helical mesh, 1st must be 0
    DATA.HELh2 = [1.0]; % For a helical mesh, last must be 1
    DATA.HELZ    = [100]; % Gives number of 'vertical' elements per h1/h2 increment per 2pi pitch
    DATA.HELdp   = 100; % Gives number of 'arc length' elements per 2pi pitch
       
    % [Z = 100, dp = 200] for conformal alpha = 15 
    
    % Hopper meshing struct
    DATA.HOPh1 =   [500  1000 2000 6000]; % Heights relative to base of hopper, only considered if hInc = 1
    DATA.HOPh2 =   [1000 2000 6000 DATA.hHtot];
    DATA.HOPZ =    [10   10   10   10];
    DATA.HOPt =    [6    6    6    5.5];
    DATA.HOPType = ['j'  'j'  'j'  'j'];
    
    % NOTE - Element types; a - S4; b - S4R; c - S4R5; d - S3; e - STRI3; f - STRI65; g - S8R; h - S8R5; j - S9R5
    % GENERAL-PURPOUSE (thick & thin) elements - S3; S4; S4R;
    % THICK elements - S8R;
    % THIN elements - STRI3; STRI65; S4R5; S8R5; S9R5; 
    % LARGE (FINITE) strain - S3; S4; S4R; 
    % SMALL strain (but finite rotations) - STRI3; S4R5; STRI65; S8R; S8R5; S9R5;   
    
    % Cylindrical meshing struct
    DATA.CYLh1   = [0  ]; 
    DATA.CYLh2   = [1  ]; 
    DATA.CYLZ    = [100]; 
    DATA.CYLt    = [t  ]; 
    DATA.CYLType = ['j']; 
    DATA.CYLRoofZ = 10;
       
    DATA.CIRCt1 = [0  ]; % Circumferential position in degrees relative to T = 0.
    DATA.CIRCt2 = [1  ]; 
    DATA.CIRCT  = [100];
        
    DATA.H_W_Toggle = 1; % Spiral welds on or off (1 - on); cylinder only - Rotter & Teng 1989 formulation
    DATA.H_W_Type = 'A'; % Only 'A' or 'B'
    DATA.H_W_Ampt = ampt; % In terms of local wall thicknesses
    
    DATA.C_W_Toggle = 0; % Circumferential welds on or off (1 - on); cylinder only
    DATA.C_W_Type = ['A'];
    DATA.C_W_Z = [0];
    DATA.C_W_Ampt = 0;
    %DATA.C_W_Type = ['A' 'A' 'A' 'A' 'A' 'A' 'A' 'A']; % Type A or B weld
    %DATA.C_W_Z = [1000 2300 3600 4800 6000 8000 10000 12000]; % Height above base of CYLINDER
    %DATA.C_W_Ampt = 0.5*1.487797589*ones(1,length(DATA.C_W_Z)); % Inward amplitudes (in local wall thicknesses) of weld depression
    %DATA.C_W_Ampt = 1*[1.487797589 1.629800601 1.629800601 1.822172467 1.822172467 2.104063529 2.104063529 2.104063529]; % Inward amplitudes (in local wall thicknesses) of weld depression

    DATA.A_W_Toggle = 0; % Axial welds on or off (1 - on); cylinder only
    DATA.A_W_Type = 'A'; % One type for all axial welds 
    DATA.A_W_Z = [0]; % Heights where welds will be displaced by the phase
    DATA.A_W_Ampt = [0.5]; % Inward amplitudes (in local wall thicknesses) of weld depression: must be 1 longer than .Z
    DATA.A_W_t = [6 5.5 5.4 5.3 5.2 5.0];
    DATA.A_W_Tfreq = 60; % Frequency around circumference (every X degrees)
    DATA.A_W_Tphase = 30; % Welds in subsequent strake out by a phase angle (in degrees)

    DATA.Super_ell = [0 3.162*3 90 2.2915562087356 1.84183448675609 11800 1]; % Super-elliptical centre-channel cylinder flattening, centered at T = 0; cylinder only
    % [Toggle (1 - on)  Amplitude of inward depression (in actual dimensional values) theta_S  power_p  power_q  Toggle sine  Z_of_max  (1 - on)]

    DATA.Circ_flat = [0 3*3.162 90 11800 1]; % Circular central flattening
    % [Toggle (1 - on)  Amplitude of inward depression (in actual dimensional values)  angular extent from channel centre (degrees)  Toggle sine  Z_of_max  (1 - on)]
    
    DATA.Axi_sin = [0 5 1]; % Inward global axial sine wave; cylinder only
    % [Toggle (1 - on)  Amplitude of inward depression (in actual dimensional values)  No. of half-waves]

    DATA.Circ_cos = [0 0.01 3 1]; % Circumferential global cosine wave; cylinder only
    % [Toggle (1 - on)  Amplitude of inward depression (in actual dimensional values)  No. of waves  Toggle sine (1 - on)]
    
    DATA.Cyl_perturbation = [0 23 2 0.01]; 
    % [Toggle (1 - on)  No. of circumferential full waves (n)  No. of axial half-waves (m)  Amplitude (in actual dimensional values)] 

    DATA.A_Type = [3]; % Any number of steps may be created as a vector, with any number of loads (assuming they're not mutually exclusive)
    % 1 = Static - General; 2 = Static - Buckle; 3 = Static - Riks
    DATA.A_GN = [1]; % Include geometric nonlinearity per step (not valid for Buckle)? (1 - yes)
    DATA.A_LBA_no = []; % ID first then no. of +ive LBA modes
    DATA.A_Riks_step = [1 stp]; % ID first then increment size
    DATA.A_Riks_inc = [1 1000]; % ID first then no. of increments
    DATA.A_Riks_max = [1 100]; % ID first then max value of load factor
    DATA.A_Restart = 0; % generate restart file for LAST step (nonlinear Riks step only)? (1 - yes)
    DATA.A_Riks_stepR = [1e-5]; % Increment size of Restart step
    DATA.A_Riks_incR = [1000]; % No. of oncrements of Restart step

    DATA.Load_Type = [1 1 10]; %...
    %2 1 10; 2 2 25; 2 3 40; ...
    %3 1 20; 3 2 45];
    % RIGID TOP only (see below)
    % 1 = point force in X; 2 = point force in Y; 3 = point force in Z ('vertical');
    % 4 = moment about X; 5 = moment about Y; 6 = moment about Z;
    % 7 = rotation about X; 8 = rotation about Y; 9 = rotation about Z;
    % RIGID BOTTOM only (see below)
    % .1 = point force in X; .2 = point force in Y; .3 = point force in Z;
    % .4 = moment about X; .5 = moment about Y; .6 = moment about Z;
    % .7 = rotation about X; .8 = rotation about Y; .9 = rotation about Z;
    % ELSE
    % 10 = Cylinder top line load
    % 15/152 = Cylinder internal uniform pressure/friction (arbitrary magnitudes required, both for cylinder only)
    % 16/162 = Hopper internal uniform pressure/friction (arbitrary magnitude required)
    % 17 = Custom circumferential/axial distributions of internal pressure (as given in the inline functions, standalone only)
    % 19 = Eccentric silo pressures, per Rotter 1986
    % 20 = Janssen concentric silo pressures; 25 = Janssen eccentric silo pressures, per EN 1991-4
    % 201 = Janssen concentric + load defined in the inline function (as in 17)
    % 30 = Reimbert concentric silo pressures; 35 = Reimbert eccentric silo pressures, per EN 1991-4 
    % NOTE - For 20, 25, 30 & 35 the arbitrary magnitudes refer to the kc value
    % 40 = Walker concentric hopper pressures
    % 45/46 = Walker eccentric hopper pressures, per EN 1991-4 and own work (45 - smooth; 46 - hypersmooth; arbitrary magnitudes required)
    % 50 = EN 1993-4-1 Appendix C wind loading on silo cylinder (isolated silo)
    % 51 = EN 1993-4-1 Appendix C wind loading on silo cylinder (silo in battery)
    % 100 = External silo pressure data from .cvv file
    % NOTE - 20, 25, 30, 35, 50, 51 & 100 in any combination in any one step is not allowed
    % NOTE - 16, 162, 40, 45 & 46 in any combination in any one step is not allowed
    DATA.Load_A_mag = [1 1 Ncr];% ...
    %2 1 1.3; ...
    %3 1 1.4];
    DATA.Load_Ears =  []; % Ears are on by default (1 - yes), choose where you DON'T want them (loads 19, 20, 201, 25, 30, 35 only) - [1 1 1] format
    DATA.Load_Factor = []; %[2 2 10; 3 2 10]; % Factor is 1 everywhere by default, choose where you want a different one
    DATA.Load_EN_Dis = []; % EN 1991-4 discharge factors are OFF by default (1 - yes), choose where you DO want them (loads 20, 201, 30 only) - [1 1 1] format
    DATA.Load_Hopper_data = []; %2 3 45 -0.5; ...
    %3 2 45 -0.5]; % For hopper eccentric discharge only [1 1 theta0 eta]; eta should be -ve for a reduction in pressures
    DATA.Load_Channel_data = [1 1 0 0.25]; % For cylinder eccentric discharge only [1 1 channel0 kc];
    DATA.Load_QVV = ['AJS_Thesis_B_silo_ECCec100']; % Define external .qvv file, no extension
    DATA.Load.Inline.F = '1';%'-0.25*cos( (Z/1e4)*(pi*Th)/(2*T0) )'; % Definition of inline test function (use Th and Z as circumferential and axial coordinates)
    DATA.Load.Inline.T0 = 15; % Theta0 (or auther auxiliary variable) for the inline function

    DATA.BC_Top   = [1 1 1 0 0 0 0];
    DATA.BC_Trans = [1 1 1 1 0 0 0];
    % 2 1 1 1 0 0 0; ...
    % 3 1 1 1 0 0 0]; % Active step; U1 (R); U2 (Theta); U3 (Z); UR1 (@R); UR2 (@Theta); UR3 (@Z); (1 is restrained) - Cylindrical CSYS
    % For symmetry to [1 0 0 1 1 1 0];
    DATA.BC_Top_RigidRing = [];
    DATA.BC_Trans_RigidRing = [];
    % Active step; U1 (X); U2 (Y); U3 (Z); UR1 (@X); UR2 (@Y); UR3 (@Z); (1 is restrained) - Cartesian CSYS
    DATA.BC_T0 = []; % Active step; 1 = x-symmetry; -1 = x-antisymmetry; 2 = y-symmetry; -2 = y-antisymmetry
    DATA.BC_Ttot = [];
    DATA.BC_Pit = []; % Hopper base boundary conditions (last bit is whether for the constraints to be inclined or not)

    DATA.Mat_Name = 'STEEL';
    DATA.Mat_Type = 1; % 1 - isotropic (E and v only); 2 - lamina (plane stress orthotropic) (E1, E2, v12, G12 and (optional) G13, G23)
    DATA.Mat_E = [E  ]; % order: E1 (E if isotropic), E2 (optional)
    DATA.Mat_v = [v];   % v (or v12 if orthotropic)
    DATA.Mat_G = [2e3]; % ignored for isotropic; for plane stress orthotropic: G12 and (optional, for transverse shear), G13, G23
    DATA.Mat_sy = [];%fy fu  fu ]; % yield stress vs. plastic strain tabular data, leave empty if no plasticity required
    DATA.Mat_ep = [];%0  eup 2*eup]; % first value of plastic strains must be zero

    DATA.Sol_Name = 'Cement'; % name of stored material for information purposes (EN 1991-4:2006 Annex E)
    DATA.Sol_Wall = 'D2'; % Wall roughness type, choose from D1, D2 or D3
    [DATA.Sol_Weight,DATA.Sol_K,DATA.Sol_Repose,DATA.Sol_Equiv,DATA.Sol_mewU,DATA.Sol_mewL,DATA.Sol_Frang,DATA.Sol_Ch,DATA.Sol_Cw] = AJSABQ_solids_database(DATA.Sol_Name,DATA.Sol_Wall,DATA.Rtot);

    DATA.Imp.Toggle = 0; % Include the imperfection key word?
    DATA.Imp.File = 'D:\SIMULIA\Jobs\AJS_CVS_LBA_var_050';
    DATA.Imp.Ampt = 5.412659; % Factoring amplitude of the imperfections
    DATA.Imp.NodeFile = 0; % Toggle whether to include the node file key word
    
    DATA.Spring.Toggle = 0; % Include elastic restraint due to solid stiffness? (1 for yes)
    DATA.Spring.Constant = 0.7; % Global spring constant (N/mm)
    DATA.Spring.Discharge = 1; % Special routine for discharge? (1 - yes)
    % NOTE; For the above, only ONE discharge load case is allowed for all steps (i.e. ONE of 19, 20, 25, 30, 35, 40, 45, 46, 100)
    % In this case, the spring stiffness will calculated as a function of the local pressure, and there will be NO springs inside the flow channel
    
    if F_save == 1; save AJS_B_Suite_Inputs_UNI_GMNIA DATA; DATA.F_save = 1; 
    else DATA.F_save = 0; 
    end
end

if F_run == 0; DATA.go = 0; disp('Calculation not requested.'); end