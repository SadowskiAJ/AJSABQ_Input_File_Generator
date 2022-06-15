% ABAQUS Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 16:57 (previously 28/03/14 - 18:17)

% Please run this generator from the command line by typing: 
% clear; MatMPI_Delete_all; eval(MPI_Run('AJSABQ_drivrutin',N,{})) 
% where N is the number of processors/threads you have/wish to assign

% Unzip MatlabMPI to be in the same directory as the .m files below.
% You must add MatlabMPI/src to the Matlab path ('selected folders and subfolders'). Make sure to save.

% The following complete list of 35 (34 compulsory, 1 optional, 4 external) files is required for successful execution:
%%%%%%% 1)  AJSABQ_drivrutin.m - DRIVING PROGRAM - DO NOT MODIFY
%%%%%%% 2)  AJSABQ_input.m - Subroutine to hold the input data
%%%%%%% 3)  AJSABQ_solids_database.m - Subroutine holding EN 1994-1 Annex E Particulate solids properties
%%%%%%% 4)  AJSABQ_structs.m - Global struct definition subroutine
%%%%%%% 5)  AJSABQ_error_check.m - Consistency checking subroutine (in progress)
%%%%%%% 6)  AJSABQ_assemble.m - General input data processing subroutine
%%%%%%% 7)  AJSABQ_qvv_process.m - External pressure .qvv file processing subroutine
%%%%%%% 8)  AJSABQ_pre_meshing.m - Preparation for meshing subroutine
%%%%%%% 9)  AJSABQ_meshing.m - Main meshing subroutine for regular shell elements
%%%%%%% 10) AJSABQ_helical_phoenix.m - Main meshing subroutine for helical shell elements
%%%%%%% 11) AJSABQ_helical_mapping - Auxiliary subroutine to the above
%%%%%%% 12) AJSABQ_triangulation.m - Auxiliary triangular meshing subroutine for helical shell elements
%%%%%%% 13) AJSABQ_Delaunay.m - Delaunay triangulation algorithm - auxiliary to the above
%%%%%%% 14) AJSABQ_springs.m - Main meshing subroutine for spring elements
%%%%%%% 15) AJSABQ_write_heading.m - Subroutine to write the Introduction to the .inp file
%%%%%%% 16) AJSABQ_write_part.m - Subroutine to write the nodes and mesh to the .inp file
%%%%%%% 17) AJSABQ_write_assembly.m - Subroutine to write node and element sets to the .inp file
%%%%%%% 18) AJSABQ_fwrite_nset.m - Auxiliary subroutine to the above
%%%%%%% 19) AJSABQ_write_material.m - Subroutine to write material definitions to the .inp file
%%%%%%% 20) AJSABQ_write_step.m - Subroutine to write the step and load definitions to the .inp file
%%%%%%% 21) AJSABQ_write_restart.m - Subroutine to write the restart step and load definitions to a separate .inp file
%%%%%%% 22) AJSABQ_fwrite_load.m - Auxiliary subroutine to the above
%%%%%%% 23) AJSABQ_objective_function.m - Auxiliary minimisation function to the above
%%%%%%% 24) AJSABQ_draw_2dplanar.m - Subroutine to plot the planar flow channel geometry (if applicable)
%%%%%%% 25) AJSABQ_draw_3dmesh.m - Subroutine to plot the structural mesh of the shell
%%%%%%% 26) AJSABQ_draw_3dsurf.m - Subroutine to plot the surface element-by-element pressures
%%%%%%% 27) AJSABQ_fvectorise.m - Auxiliary subroutine to the above
%%%%%%% 28) AJSABQ_summary.m - Subroutine to display a summary of the model to the Matlab workspace window
%%%%%%% 29) AJSABQ_paths.m - Subroutine to specify the paths for ABAQUS CAE data extraction
%%%%%%% 30) AJSABQ_fwrite_python.m - Auxiliary subroutine to the above: writes a Python script for ABAQUS CAE data extraction
%%%%%%% 31) EXTERNAL_unique_no_sort.m - External subroutine to remove duplicate vector entires without sorting the vector (C) 2007, Caitlin Bever
%%%%%%% 32) EXTERNAL_inpoly.m - External 'point-in-polygon' subroutine (C) 2007, Darren Engwirda
%%%%%%% 33) EXTERNAL_nearest_neighbour.m - External 'nearest-neighbour' subroutine - (C) 2006, Richard Brown
%%%%%%% 34) EXTERNAL_notify.wav - External sound file to signal the successful generation of the input file - (C) Microsoft
% If the eccentric discharge silo pressures are required, to make the .qvv file, an independent .m file generates them:
%%%%%%% 35) AJSABQ_silo_pressures.m - External program to generate the pressure .qvv file

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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hereafter please do not modify the code unless you know exactly what you are doing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

PP = struct('CPUs',[],'comm',[]); % Parallel execution struct
TEMP = struct('ddd1',[],'ddt1',[],'go',[],'fid',[],'fidR',[]); % Temporary data struct
clc; MPI_Init; PP.comm = MPI_COMM_WORLD; PP.CPUs = MPI_Comm_size(PP.comm); PP.rank = MPI_Comm_rank(PP.comm); % Initialising MPI
PP.master = 0; PP.slave = [1:(PP.comm.size)-1]; if PP.CPUs == 1; PP.slave = 999; end
TEMP.ddd1 = datestr(now); TEMP.ddt1 = clock; 

if isempty(find(PP.master == PP.rank,1)) == 0 % If master thread  
    [GEO,NAME,DRAW,MATERIAL,RSL,STEP,SOLID,PATHS,SPRING] = AJSABQ_structs(); % Loading global structs necessary for initial data processing

    [DATA] = AJSABQ_input(); % Loading model input data from separate function
    TEMP.go = DATA.go;
 
    if TEMP.go ~= 0
        GEO.cyl.Inc = DATA.cInc; GEO.hop.Inc = DATA.hInc; GEO.roof.Inc = DATA.rInc;
        NAME.FileName = DATA.Name; NAME.Heading = DATA.Heading; NAME.Description = DATA.Description;
        GEO.Rtot = DATA.Rtot; GEO.cyl.Htot = DATA.cHtot; GEO.hop.Htot = DATA.hHtot; GEO.roof.Angle = DATA.rAngle; GEO.roof.Thick = DATA.rThick;
        RSL.HOP.h1 = DATA.HOPh1; RSL.HOP.h2 = DATA.HOPh2; RSL.HOP.Z = DATA.HOPZ; RSL.HOP.t = DATA.HOPt; RSL.HOP.Type = DATA.HOPType;
        RSL.CYL.h1 = DATA.CYLh1; RSL.CYL.h2 = DATA.CYLh2; RSL.CYL.Z = DATA.CYLZ; RSL.CYL.t = DATA.CYLt; RSL.CYL.Type = DATA.CYLType; RSL.CYL.RoofZ = DATA.CYLRoofZ;
        RSL.HEL.h1 = DATA.HELh1; RSL.HEL.h2 = DATA.HELh2; RSL.HEL.Z = DATA.HELZ; RSL.HEL.dp = DATA.HELdp;
        RSL.CIRC.t1 = DATA.CIRCt1; RSL.CIRC.t2 = DATA.CIRCt2; RSL.CIRC.T  = DATA.CIRCT;
        GEO.feature.H_Weld.Toggle = DATA.H_W_Toggle; GEO.feature.H_Weld.Type = DATA.H_W_Type; GEO.feature.H_Weld.Ampt = DATA.H_W_Ampt;
        GEO.feature.C_Weld.Toggle = DATA.C_W_Toggle; GEO.feature.C_Weld.Type = DATA.C_W_Type; GEO.feature.C_Weld.Z = DATA.C_W_Z; GEO.feature.C_Weld.Ampt = DATA.C_W_Ampt;
        GEO.feature.A_Weld.Toggle = DATA.A_W_Toggle; GEO.feature.A_Weld.Type = DATA.A_W_Type; GEO.feature.A_Weld.Z = DATA.A_W_Z; GEO.feature.A_Weld.Ampt = DATA.A_W_Ampt;
        GEO.feature.A_Weld.t = DATA.A_W_t; GEO.feature.A_Weld.Tfreq = DATA.A_W_Tfreq; GEO.feature.A_Weld.Tphase = DATA.A_W_Tphase;
        GEO.feature.Super_ell = DATA.Super_ell; GEO.feature.Circ_flat = DATA.Circ_flat;
        GEO.feature.Axi_sin = DATA.Axi_sin; GEO.feature.Circ_cos = DATA.Circ_cos;
        GEO.feature.Cyl_perturbation = DATA.Cyl_perturbation;
        GEO.helix.Toggle = DATA.helical.Toggle; GEO.helix.choice = DATA.helical.choice;
        GEO.helix.pitch = DATA.helical.pitch; GEO.helix.incline = DATA.helical.incline; 
        GEO.helix.element = DATA.helical.element; GEO.helix.thick = DATA.helical.thick;
        GEO.helix.conformal = DATA.helical.conformal;
        DRAW.surf_3D.step = DATA.surf_3D.step; DRAW.surf_3D.which = DATA.surf_3D.which; DRAW.planar_2D.step = DATA.planar_2D.step;
        PATHS.Toggle = DATA.paths.Toggle; PATHS.Python.Toggle = DATA.python.Toggle; PATHS.Python.Copy = DATA.python.Copy;
        DRAW.mesh_3D.toggle = DATA.mesh_3D.toggle; DRAW.mesh_3D.nodes = DATA.mesh_3D.nodes; DRAW.mesh_3D.elements = DATA.mesh_3D.elements;
        STEP.Ana.Type = DATA.A_Type; STEP.Ana.GN = DATA.A_GN; STEP.Ana.LBA_no = DATA.A_LBA_no; STEP.Ana.Riks_step = DATA.A_Riks_step; STEP.Ana.Riks_max = DATA.A_Riks_max;
        STEP.Ana.Riks_inc = DATA.A_Riks_inc; STEP.Ana.Restart = DATA.A_Restart; STEP.Ana.Riks_stepR = DATA.A_Riks_stepR; STEP.Ana.Riks_incR = DATA.A_Riks_incR;
        STEP.Load.Type = DATA.Load_Type; STEP.Load.A_mag = DATA.Load_A_mag; STEP.Load.Ears = DATA.Load_Ears; STEP.Load.Factor = DATA.Load_Factor; STEP.Load.EN_Dis = DATA.Load_EN_Dis;
        STEP.Load.Hopper_data = DATA.Load_Hopper_data; STEP.Load.Channel_data = DATA.Load_Channel_data; STEP.Load.QVV = DATA.Load_QVV;
        STEP.BC.Top = DATA.BC_Top; STEP.BC.Top_RigidRing = DATA.BC_Top_RigidRing; STEP.BC.Trans = DATA.BC_Trans; STEP.BC.Trans_RigidRing = DATA.BC_Trans_RigidRing; 
        STEP.BC.T0 = DATA.BC_T0; STEP.BC.Ttot = DATA.BC_Ttot; STEP.BC.Pit = DATA.BC_Pit;
        MATERIAL.Name = DATA.Mat_Name; MATERIAL.Type = DATA.Mat_Type; MATERIAL.E = DATA.Mat_E; MATERIAL.v = DATA.Mat_v; MATERIAL.G = DATA.Mat_G; MATERIAL.sy = DATA.Mat_sy; MATERIAL.ep = DATA.Mat_ep;
        SOLID.Name = DATA.Sol_Name; SOLID.Wall = DATA.Sol_Wall; SOLID.Weight = DATA.Sol_Weight; SOLID.K = DATA.Sol_K; SOLID.Repose = DATA.Sol_Repose; SOLID.Equiv = DATA.Sol_Equiv;
        SOLID.mewU = DATA.Sol_mewU; SOLID.mewL = DATA.Sol_mewL; SOLID.Ch = DATA.Sol_Ch; SOLID.Cw = DATA.Sol_Cw; SOLID.Frang = DATA.Sol_Frang;
        SPRING.toggle = DATA.Spring.Toggle; SPRING.constant = DATA.Spring.Constant; SPRING.Discharge = DATA.Spring.Discharge;

        [TEMP.go] = AJSABQ_error_check(GEO,NAME,DRAW,MATERIAL,RSL,STEP,DATA,SOLID,SPRING,PATHS,PP); % Subroutine to check the consistency of the input data
    end
 
    if TEMP.go ~= 1
        disp(' '); disp(['      Generation failed. ']);
    else
        OP = struct('R',[],'H',[],'Hc',[],'Hh',[],'Htot',[],'T',[],'Hpitch',[],'N_p',[],'cvv',[],'cvv_file',[],'Sol',[],'Con',[],'ConH',[],... % Operational struct
            'i1',[],'i2',[],'e1',[],'e2',[],'Hmin',[],'Hmax',[],'div',[],'RemTOT',[],'Hod',[],'Class',[],...
            'type',[],'mag',[],'hop_th0',[],'hop_eta',[],'ears',[],'factor',[],'EN_Dis',[],'kc',[],'chan0',[],'Riks',[],...
            'count_1',[],'count_2',[],'count_3',[],'plastic',[],'Hop',[],'Cyl',[],'Roof',[],'CWeld',[],'Inline',[],'D1',[],'D2',[],'xoff',[],...
            'Mat_Save',[],'Copy',[],'Imp',[],'Heltoggle',[],'Helt',[],'Heltype',[],'Helorder',[],'Hel_Quad',[],'Hel_Tri',[],'Tols',[]);
        OP.Con = struct('Ar',[],'U',[],'zo',[],'pho',[],'nR',[]);
        OP.ConH = struct('F',[],'n',[],'pvft',[]);
        OP.Ecc = struct('ec',[],'THc',[],'Psi',[],'Uwc',[],'Usc',[],'Ac',[]);
        OP.Inline = struct('S',[],'F',[],'T0',[]);
        OP.Inline.S = DATA.Load.Inline.F; OP.Inline.F = inline(DATA.Load.Inline.F,'Th','T0','Z'); OP.Inline.T0 = DATA.Load.Inline.T0; OP.Mat_Save = DATA.F_save;
        OP.Copy = struct('toggle',[],'destination',[]);
        OP.Copy.toggle = DATA.copy.Toggle; OP.Copy.destination = DATA.copy.Destination;
        OP.Imp = struct('toggle',[],'file',[],'ampt',[],'nodefile',[]);
        OP.Imp.toggle = DATA.Imp.Toggle; OP.Imp.file = DATA.Imp.File; OP.Imp.ampt = DATA.Imp.Ampt; OP.Imp.nodefile = DATA.Imp.NodeFile; OP.Steps = size(DATA.A_Type,2);

        [GEO,RSL,STEP,OP,NAME] = AJSABQ_assemble(GEO,RSL,STEP,SOLID,OP,NAME,MATERIAL); % Subroutine to process the input data

        QVV = struct('Heading',[],'Z',[],'T',[],'qn',[],'SOLID',[],'CHAN',[],'dTT',[],'File',[]); % Struct to hold data from .qvv file
        QVV.SOLID = struct('gam',[],'mu_i',[],'mu_w',[],'Ks_w',[],'Ks_c',[],'Kc_w',[],'Kc_s',[],'factor',[]);
        QVV.CHAN = struct('ec',[],'hJ',[],'h0',[],'power',[]);
        
        if OP.cvv == 1 % If external .qvv file is requested
            [QVV,TEMP.go] = AJSABQ_qvv_process(QVV,OP); NAME.Heading = QVV.Heading; QVV.File = OP.cvv_file; % Subroutine to process the external .qvv file
        end

        if GEO.helix.Toggle ~= 1 % Rectangular meshing structs
            NODES = struct('index',[],'R',[],'T',[],'Z',[]); % Struct listing node index and radius, theta and z coordinates (cylindrical polar system)
            ELEMENTS = struct('index',[],'nodes',[],'region',[],'type',[],'t',[],'alt',[]); % Struct listing shell element indices and node coordinates plus other information
            NSET = struct('Bottom',[],'Top',[],'Right',[],'Left',[],'Czelusc',[],'Riks',[],'Axis',[],'Path',[]); % Struct listing node sets
            NSET.Path = struct('Centre',[],'Edge',[],'Opposite',[]); % Struct listing node paths
            ELSET = [];
        elseif GEO.helix.Toggle == 1 % Helical meshing structs
            NODES = struct('index',[],'x',[],'y',[],'X',[],'Y',[],'Z',[]);
            ELEMENTS = struct('index',[],'nodes',[],'type',[]); NSET = [];
            OP.Tols.tolBG = 1e-5; % 'Large' tolerance (1e-5)
            OP.Tols.tolNN = 1e-10; % Nearest neighbour tolerance (1e-10)
            OP.Tols.tolRND = 20; % Round to this no. of decimal points (20)
            OP.Tols.tolH = 1e-10; % Vertical tolerance (1e-10)
            OP.Tols.tolS = 1e-10; % Side tolerance (1e-10)
            OP.Tols.sideE = 1; % Side element gap (2)
            OP.Tols.interior = 0; % Include interior nodes in side regions?
            OP.Tols.scale = 1; % Scale side regions to [0,1] 2D domain?
        end
        
        [OP] = AJSABQ_pre_meshing(RSL,OP,PP); % Pre-meshing subroutine
    end
    MPI_Bcast(PP.master,3999,PP.comm,TEMP.go);
end

if isempty(find(PP.slave == PP.rank,1)) == 0; [TEMP.go] = MPI_Recv(PP.master,3999,PP.comm); end
if TEMP.go == 0 && isempty(find(PP.slave == PP.rank,1)) == 0; exit; end
if TEMP.go == 1 
    if isempty(find(PP.master == PP.rank,1)) == 0; MPI_Bcast(PP.master,4000,PP.comm,GEO,RSL,OP,MATERIAL,STEP,NAME,NODES,ELEMENTS,NSET,SOLID,QVV,DRAW,SPRING); end
    if isempty(find(PP.slave == PP.rank,1)) == 0; [GEO, RSL, OP, MATERIAL, STEP, NAME, NODES, ELEMENTS, NSET, SOLID, QVV, DRAW, SPRING] = MPI_Recv(PP.master,4000,PP.comm); end
    if (isempty(find(PP.master == PP.rank,1)) == 0 || isempty(find(PP.slave == PP.rank,1)) == 0) && TEMP.go == 1
        if GEO.helix.Toggle ~= 1
            % Meshing regular shell elements
            [NODES,ELEMENTS,NSET,OP] = AJSABQ_meshing(PP,RSL,OP,NODES,ELEMENTS,NSET,GEO,MATERIAL);
        end
        if GEO.helix.Toggle == 1
            % Meshing helical shell elements
            [NODES,ELEMENTS,OP] = AJSABQ_helical_phoenix(PP,RSL,OP,NODES,ELEMENTS,GEO,MATERIAL);
            [NODES,ELEMENTS,OP,NSET,ELSET] = AJSABQ_triangulation(PP,RSL,OP,NODES,ELEMENTS,GEO,MATERIAL);
        end
    end
    if isempty(find(PP.master == PP.rank,1)) == 0 && isempty(OP.RemTOT) == 0; OP.RemTOT(length(OP.RemTOT)) = []; NODES.index(OP.RemTOT) = []; NODES.R(OP.RemTOT) = []; NODES.Z(OP.RemTOT) = []; NODES.T(OP.RemTOT) = []; end

    if PP.comm.size > 2 % Only 2 processors are necessary now
        PP.comm.size = 2; PP.comm.group = [0 1]; PP.comm.machine_id = [1 1]; PP.comm.machine_db.n_proc = 2; PP.comm.machine_db.id_stop = 2; PP.slave = 1;
        if PP.rank > 1; exit; end
    end
    if GEO.helix.Toggle == 1
        if PP.rank > 0; exit; end
    end
    
    % Writing the ABAQUS input file
    % Opening
    if PP.rank == PP.master
        clc; disp(' '); disp('      Writing ABAQUS input file'); disp(' '); disp(['      Number of assigned processors/threads:    ',num2str(PP.CPUs)]);
        TEMP.fid = fopen(NAME.FileName,'w');
    end
   
    if (STEP.Ana.Restart == 1 && PP.rank == PP.slave) || (STEP.Ana.Restart == 1 && PP.CPUs == 1)
        TEMP.fidR = fopen(NAME.RestartName,'w'); 
    end

    if isempty(find(PP.master == PP.rank,1)) == 0
        
        % Computing paths
        if PATHS.Toggle == 1
            [PATHS] = AJSABQ_paths(PATHS,NODES,GEO,OP,NAME,RSL.Ttot);
        end

        % Writing Introduction
        [TEMP.fid] = AJSABQ_write_heading(TEMP.fid,STEP,NAME,GEO,PP,OP,SOLID,RSL,NODES,ELEMENTS.index,PATHS,SPRING);

        % Writing Part
        [TEMP.fid] = AJSABQ_write_part(TEMP.fid,RSL,NODES,ELEMENTS,MATERIAL,OP);

        % Writing Assembly
        [TEMP.fid,NSET] = AJSABQ_write_assembly(TEMP.fid,GEO,OP,RSL,STEP,NSET,ELSET,NODES,ELEMENTS,SPRING);
        
        % Writing Material
        [TEMP.fid] = AJSABQ_write_material(TEMP.fid,MATERIAL);
        
        % Writing Step & Loads
        for S = 1:STEP.Tot 
            ST_P = struct('R',[],'T',[],'Z',[]); ST_F = struct('R',[],'T',[],'Z',[]);
            pos = find(DRAW.surf_3D.step == S, 1);
            if ~isempty(pos); DRAW3DS = DRAW.surf_3D.which(S); else DRAW3DS = 0; end
            pos2 = find(DRAW.planar_2D.step == S, 1);
            if ~isempty(pos2); DRAW2DP = 1; else DRAW2DP = 0; end
            [TEMP.fid,OP,~,~] = AJSABQ_write_step(TEMP.fid,S,STEP,OP,SOLID,ELEMENTS,NODES,QVV,DRAW3DS,DRAW2DP,ST_P,ST_F,SPRING.toggle); 
            clear('ST_P','ST_F','pos','pos2','DRAW2DP','DRAW3DS');
        end
        clear('S');
        
        % Plotting mesh
        if DRAW.mesh_3D.toggle == 1; AJSABQ_draw_3dmesh(ELEMENTS,NODES,SPRING,OP,DRAW.mesh_3D,RSL.Ttot); end     
    end
    
    % Writing Restart File
    if PP.rank == PP.master && STEP.Ana.Restart == 1 && PP.CPUs ~= 1; MPI_Send(PP.slave,40000,PP.comm,NODES,ELEMENTS); fclose(TEMP.fid); end
    if PP.rank == PP.slave && STEP.Ana.Restart == 1 && PP.CPUs ~= 1; [NODES, ELEMENTS] = MPI_Recv(PP.master,40000,PP.comm); end
    if (STEP.Ana.Restart == 1 && PP.rank == PP.slave) || (STEP.Ana.Restart == 1 && PP.CPUs == 1)
        ST_P = struct('R',[],'T',[],'Z',[]); ST_F = struct('R',[],'T',[],'Z',[]);
        pos = find(DRAW.surf_3D.step == STEP.Tot);
        if ~isempty(pos); DRAW3DS = DRAW.surf_3D.which(STEP.Tot); else DRAW3DS = 0; end
        pos2 = find(DRAW.planar_2D.step == STEP.Tot);
        if ~isempty(pos2); DRAW2DP = 1; else DRAW2DP = 0; end
        [TEMP.fidR] = AJSABQ_write_restart(TEMP.fidR,STEP.Tot,STEP,OP,SOLID,ELEMENTS,NODES,QVV,DRAW3DS,DRAW2DP,ST_P,ST_F,SPRING.toggle,PP.CPUs);
        fclose(TEMP.fidR); MPI_Send(PP.master,40001,PP.comm,1);
    end
    
    % Terminating & displaying summary to Matlab window
    if PP.rank == PP.master
        if (PP.CPUs == 1 && STEP.Ana.Restart == 1) || STEP.Ana.Restart ~= 1; fclose(TEMP.fid); end
        if STEP.Ana.Restart == 1 && PP.CPUs ~= 1; a0 = MPI_Recv(PP.slave,40001,PP.comm); 
            if a0 == 1; AJSABQ_summary(TEMP,STEP,RSL,NAME,OP,DRAW,length(NODES.index),length(ELEMENTS.index),GEO,PATHS,SPRING.toggle,length(SPRING.index),SOLID); clear('a0'); end
        else
            AJSABQ_summary(TEMP,STEP,RSL,NAME,OP,DRAW,length(NODES.index),length(ELEMENTS.index),GEO,PATHS,SPRING.toggle,length(SPRING.index),SOLID);
        end
    end
    if PP.rank == PP.slave; exit; end
    
    % Final Clean-up
    if OP.Copy.toggle == 1 
        [a1,b1,c1]=copyfile(NAME.FileName,OP.Copy.destination,'f');
        if a1 == 0; disp('      WARNING: Main .inp file not copied.'); 
        else disp(['      Main file copied to ',OP.Copy.destination]); clear('a1','b1','c1');
        end
        if ~isempty(NAME.RestartName); [a2,b2,c2]=copyfile(NAME.RestartName,OP.Copy.destination,'f');
            if a2 == 0; disp('      WARNING: Restart .inp file not copied.'); 
            else disp(['      Restart file copied to ',OP.Copy.destination]); clear('a2','b2','c2');
            end
        end
    end
    if PATHS.Python.Copy == 1
        if PATHS.Python.Result == 0; disp('      WARNING: Python .py script not copied.'); 
        else disp(' '); disp(['      Python script copied to ',OP.Copy.destination]);
        end
    end
end
disp(' '); clear('ans'); disp('      Completed.');
% Thank you.