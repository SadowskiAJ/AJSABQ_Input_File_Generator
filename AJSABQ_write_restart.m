function [fidR] = AJSABQ_write_restart(fidR,S,STEP,OP,SOLID,ELEMENTS,NODES,QVV,DRAW3DS,DRAW2DP,ST_P,ST_F,SPTOG,CPUs)
% Function to write the Restart step and Load definitions of the .inp file
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:40 (previously 09/08/11 - 18:48)

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

% This step must be a nonlinear Riks step
disp(['Writing restart step (static - riks)']);

fprintf(fidR,'%s\n',['** Restart .inp file created on ',datestr(now),' with the ABQAJS input file generator Version 2.0 using ',num2str(CPUs),' processors']);
fprintf(fidR,'%s\n','**');
fprintf(fidR,'%s\n',['*RESTART, READ, STEP=',num2str(S),', INC=X, ENDSTEP']);
fprintf(fidR,'%s\n','**');
fprintf(fidR,'%s\n',['*STEP, NAME=RESTART_STATIC_NONLINEAR_RIKS, INC=',num2str(STEP.Ana.Riks_incR)]);
fprintf(fidR,'%s\n','**');
fprintf(fidR,'%s\n','*STATIC, RIKS');
fprintf(fidR,'%s\n',[num2str(STEP.Ana.Riks_stepR),',1,1e-015,',num2str(STEP.Ana.Riks_stepR),',1e5,RIKS,3,1e5']);
fprintf(fidR,'%s\n','**');

fprintf(fidR,'%s\n','*BOUNDARY, OP=NEW'); 
if isempty(STEP.BC.Top) ~= 1; for C = 2:7; if STEP.BC.Top(S,C) == 1; fprintf(fidR,'%s\n',['TOP, ',num2str(C-1),', ',num2str(C-1)]); end; end; end
if isempty(STEP.BC.Trans) ~= 1; for C = 2:7; if STEP.BC.Trans(S,C) == 1; fprintf(fidR,'%s\n',['TRANSITION, ',num2str(C-1),', ',num2str(C-1)]); end; end; end
if isempty(STEP.BC.Pit) ~= 1; for C = 2:7; if STEP.BC.Pit(S,C) == 1; fprintf(fidR,'%s\n',['PIT, ',num2str(C-1),', ',num2str(C-1)]); end; end; end
if STEP.BC.Ttot(S,2) == 1; fprintf(fidR,'%s\n','LEFT, XSYMM');
elseif STEP.BC.Ttot(S,2) == -1; fprintf(fidR,'%s\n','LEFT, XASYMM');
elseif STEP.BC.Ttot(S,2) == 2; fprintf(fidR,'%s\n','LEFT, YSYMM');
elseif STEP.BC.Ttot(S,2) == -2; fprintf(fidR,'%s\n','LEFT, YASYMM');
end
if STEP.BC.T0(S,2) == 1; fprintf(fidR,'%s\n','RIGHT, XSYMM');
elseif STEP.BC.T0(S,2) == -1; fprintf(fidR,'%s\n','RIGHT, XASYMM');
elseif STEP.BC.T0(S,2) == 2; fprintf(fidR,'%s\n','RIGHT, YSYMM');
elseif STEP.BC.T0(S,2) == -2; fprintf(fidR,'%s\n','RIGHT, YASYMM');
end
if SPTOG == 1; fprintf(fidR,'%s\n','AXIS, PINNED'); end
fprintf(fidR,'%s\n','**');

SP = 0; SF = 0;
for T = 1:STEP.Sub(S) % Writing Loads
    [fidR,OP,ST_P,ST_F,SP,SF] = AJSABQ_fwrite_load(fidR,S,T,OP,SOLID,ELEMENTS,NODES,QVV,DRAW3DS,DRAW2DP,ST_P,ST_F,SP,SF);
end

fprintf(fidR,'%s\n','**');
fprintf(fidR,'%s\n','*CONTROLS, PARAMETERS=TIME INCREMENTATION');
fprintf(fidR,'%s\n',' , , , , , , ,20, , '); fprintf(fidR,'%s\n','**');
fprintf(fidR,'%s\n','*OUTPUT, FIELD'); fprintf(fidR,'%s\n','*NODE OUTPUT');
fprintf(fidR,'%s\n','RF, TF, U, COORD');
fprintf(fidR,'%s\n','*ELEMENT OUTPUT, DIRECTIONS=YES');
fprintf(fidR,'%s\n',['S, SSAVG, E, PE, PEEQ, EE, IE, PEMAG, NE, LE, SE']);
fprintf(fidR,'%s\n','*OUTPUT, HISTORY');
fprintf(fidR,'%s\n','*NODE OUTPUT, NSET=RIKS'); fprintf(fidR,'%s\n','U1, U2, U3');
fprintf(fidR,'%s\n','**');
fprintf(fidR,'%s\n','*END STEP'); 