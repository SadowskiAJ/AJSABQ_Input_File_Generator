function [fid,OP,ST_P,ST_F] = AJSABQ_write_step(fid,S,STEP,OP,SOLID,ELEMENTS,NODES,QVV,DRAW3DS,DRAW2DP,ST_P,ST_F,SPTOG)
% Function to write the Step and Load definitions of the .inp file 
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:40 (previously 01/07/12 - 20:47)

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

fprintf(fid,'%s\n','**'); s = STEP.Ana.Type(S);
if s == 1 % STATIC - GENERAL STEP
    disp(['Writing step ',num2str(S),'/',num2str(STEP.Tot),' (static - general)']);
    OP.count_1 = OP.count_1 + 1;
    if STEP.Ana.GN(S) == 1
        for I = 1:size(STEP.Ana.Riks_step,1); if STEP.Ana.Riks_step(I,1) == S; stp = STEP.Ana.Riks_step(I,2); end; end
        non = 'NON'; nlg = ', NLGEOM=YES'; vec = [stp, 1.0, 1e-15, stp];
    else
        non = ''; nlg = ', NLGEOM=NO'; vec = [0.5, 1.0, 1e-15, 0.5];
    end
    fprintf(fid,'%s\n',['*STEP, NAME=STATIC_',non,'LINEAR_',num2str(OP.count_1),nlg]);
    fprintf(fid,'%s\n','*STATIC');
    fprintf(fid,'%1.3f, %1.f, %1.e, %1.3f\n',vec);
elseif s == 2 % STATIC - PERTURBATION STEP
    disp(['Writing step ',num2str(S),'/',num2str(STEP.Tot),' (static - perturbation)']);
    OP.count_2 = OP.count_2 + 1;
    fprintf(fid,'%s\n',['*STEP, NAME=STATIC_BUCKLE_',num2str(OP.count_2),', PERTURBATION']);
    fprintf(fid,'%s\n','*BUCKLE, EIGENSOLVER=LANCZOS');
    for I = 1:size(STEP.Ana.LBA_no,1); if STEP.Ana.LBA_no(I,1) == S; nos = STEP.Ana.LBA_no(I,2); end; end
    fprintf(fid,'%s\n',[num2str(nos),', 0, , , ,']);
elseif s == 3 % STATIC - RIKS STEP
    disp(['Writing step ',num2str(S),'/',num2str(STEP.Tot),' (static - riks)']);
    OP.count_3 = OP.count_3 + 1;
    if STEP.Ana.GN(S) == 1; non = 'NON'; nlg = ', NLGEOM=YES, '; vec = []; else non = ''; nlg = ', NLGEOM=NO, '; end
    for I = 1:size(STEP.Ana.Riks_inc,1); if STEP.Ana.Riks_inc(I,1) == S; inc = STEP.Ana.Riks_inc(I,2); end; end
    fprintf(fid,'%s\n',['*STEP, NAME=STATIC_',non,'LINEAR_RIKS_',num2str(OP.count_3),nlg,'INC=',num2str(inc)]);
    fprintf(fid,'%s\n','*STATIC, RIKS');
    for I = 1:size(STEP.Ana.Riks_step,1); if STEP.Ana.Riks_step(I,1) == S; step = STEP.Ana.Riks_step(I,2); end; end
    for I = 1:size(STEP.Ana.Riks_max,1); if STEP.Ana.Riks_max(I,1) == S; LP_max = STEP.Ana.Riks_max(I,2); end; end
    fprintf(fid,'%s\n',[num2str(step),',1,1e-015,',num2str(step),',',num2str(LP_max),',RIKS,3,1e5']);
end
fprintf(fid,'%s\n','**');

fprintf(fid,'%s\n','*BOUNDARY, OP=NEW');
if isempty(STEP.BC.Top) ~= 1; for C = 2:7; if STEP.BC.Top(S,C) == 1; fprintf(fid,'%s\n',['TOP, ',num2str(C-1),', ',num2str(C-1)]); end; end; end
if isempty(STEP.BC.Top_RigidRing) ~= 1; for C = 2:7; if STEP.BC.Top_RigidRing(S,C) == 1; fprintf(fid,'%s\n',['RIGID_NODE_TOP, ',num2str(C-1),', ',num2str(C-1)]); end; end; end
if isempty(STEP.BC.Trans_RigidRing) ~= 1; for C = 2:7; if STEP.BC.Trans_RigidRing(S,C) == 1; fprintf(fid,'%s\n',['RIGID_NODE_BOTTOM, ',num2str(C-1),', ',num2str(C-1)]); end; end; end
if isempty(STEP.BC.Trans) ~= 1; for C = 2:7; if STEP.BC.Trans(S,C) == 1; fprintf(fid,'%s\n',['TRANSITION, ',num2str(C-1),', ',num2str(C-1)]); end; end; end
if isempty(STEP.BC.Pit) ~= 1; for C = 2:7; if STEP.BC.Pit(S,C) == 1; fprintf(fid,'%s\n',['PIT, ',num2str(C-1),', ',num2str(C-1)]); end; end; end
if OP.Heltoggle ~= 1
    if STEP.BC.Ttot(S,2) == 1; fprintf(fid,'%s\n','LEFT, XSYMM');
    elseif STEP.BC.Ttot(S,2) == -1; fprintf(fid,'%s\n','LEFT, XASYMM');
    elseif STEP.BC.Ttot(S,2) == 2; fprintf(fid,'%s\n','LEFT, YSYMM');
    elseif STEP.BC.Ttot(S,2) == -2; fprintf(fid,'%s\n','LEFT, YASYMM');
    end
    if STEP.BC.T0(S,2) == 1; fprintf(fid,'%s\n','RIGHT, XSYMM');
    elseif STEP.BC.T0(S,2) == -1; fprintf(fid,'%s\n','RIGHT, XASYMM');
    elseif STEP.BC.T0(S,2) == 2; fprintf(fid,'%s\n','RIGHT, YSYMM');
    elseif STEP.BC.T0(S,2) == -2; fprintf(fid,'%s\n','RIGHT, YASYMM');
    end
end
if SPTOG == 1; fprintf(fid,'%s\n','AXIS, PINNED'); end
fprintf(fid,'%s\n','**');

SP = 0; SF = 0;
for T = 1:STEP.Sub(S) % Writing Loads
    [fid,OP,ST_P,ST_F,SP,SF] = AJSABQ_fwrite_load(fid,S,T,OP,SOLID,ELEMENTS,NODES,QVV,DRAW3DS,DRAW2DP,ST_P,ST_F,SP,SF);
end

fprintf(fid,'%s\n','**');
if s == 1 % STATIC - GENERAL STEP
    fprintf(fid,'%s\n','*OUTPUT, FIELD'); fprintf(fid,'%s\n','*NODE OUTPUT');
    fprintf(fid,'%s\n','RF, TF, U, COORD');
    fprintf(fid,'%s\n','*ELEMENT OUTPUT, DIRECTIONS=YES');
    if OP.Heltoggle == 1 && OP.Helorder == 2;
        fprintf(fid,'%s\n','S, SSAVG, E, PE, PEEQ, EE, IE, PEMAG, SE, STH');
    else
        fprintf(fid,'%s\n','S, SSAVG, E, PE, PEEQ, EE, IE, PEMAG, NE, LE, SE, STH');
    end 
elseif s == 2 % STATIC - PERTURBATION STEP
    fprintf(fid,'%s\n','*NODE FILE, GLOBAL=YES, LAST MODE=1');
    fprintf(fid,'%s\n','U');
elseif s == 3 % STATIC - RIKS STEP
    if STEP.Ana.Restart == 1; fprintf(fid,'%s\n','*RESTART, WRITE, FREQUENCY=1'); fprintf(fid,'%s\n','**'); end
    fprintf(fid,'%s\n','*CONTROLS, PARAMETERS=TIME INCREMENTATION');
    fprintf(fid,'%s\n',' , , , , , , ,20, , '); fprintf(fid,'%s\n','**');
    fprintf(fid,'%s\n','*OUTPUT, FIELD'); fprintf(fid,'%s\n','*NODE OUTPUT');
    fprintf(fid,'%s\n','RF, TF, U, COORD');
    fprintf(fid,'%s\n','*ELEMENT OUTPUT, DIRECTIONS=YES');
    if OP.Heltoggle == 1 && OP.Helorder == 2
        fprintf(fid,'%s\n','S, SSAVG, E, PE, PEEQ, EE, IE, PEMAG, SE, STH');
    else
        fprintf(fid,'%s\n','S, SSAVG, E, PE, PEEQ, EE, IE, PEMAG, NE, LE, SE, STH');
    end 
    fprintf(fid,'%s\n','*OUTPUT, HISTORY');
    if sum(STEP.BC.Top_RigidRing) > 1
        fprintf(fid,'%s\n','*NODE OUTPUT, NSET=RIGID_NODE_TOP'); fprintf(fid,'%s\n','U1, U2, U3, UR1, UR2, UR3');
    elseif sum(STEP.BC.Trans_RigidRing) > 1
        fprintf(fid,'%s\n','*NODE OUTPUT, NSET=RIGID_NODE_BOTTOM'); fprintf(fid,'%s\n','U1, U2, U3, UR1, UR2, UR3');
    else
        fprintf(fid,'%s\n','*NODE OUTPUT, NSET=RIKS'); fprintf(fid,'%s\n','U1, U2, U3');
    end
end

if OP.Imp.nodefile == 1 && s ~= 2
    fprintf(fid,'%s\n','**'); fprintf(fid,'%s\n','*NODE FILE, GLOBAL=YES, LAST MODE=1'); fprintf(fid,'%s\n','U');
end
fprintf(fid,'%s\n','**');
fprintf(fid,'%s\n','*END STEP');
if S < STEP.Tot; fprintf(fid,'%s\n','**'); end

if DRAW3DS ~= 0; AJSABQ_draw_3dsurf(ST_P,ST_F,OP,DRAW3DS,S); end