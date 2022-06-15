function [PATHS] = AJSABQ_fwrite_python(PATHS,OPCopy,NAME)
% Function to write paths to a Python script for use in ABAQUS CAE v6.9.1
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:24 (previously 24/11/09 - 16:55)

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

disp('Writing Python script.');
FileO = NAME.FileName; FileO(length(FileO)-3:length(FileO)) = []; FileO(length(FileO)+1:length(FileO)+4) = '.odb';
FileN = NAME.FileName; FileN(length(FileN)-3:length(FileN)) = []; FileN(length(FileN)+1:length(FileN)+3) = '.py';
FileT_A = NAME.FileName;
FileT_A(length(FileT_A)-3:length(FileT_A)) = [];
FileT_A(length(FileT_A)+1:length(FileT_A)+10) = '_AXIAL.tx2';
FileT_C = NAME.FileName;
FileT_C(length(FileT_C)-3:length(FileT_C)) = [];
FileT_C(length(FileT_C)+1:length(FileT_C)+20) = '_CIRCUMFERENTIAL.tx2';

% Writing file heading
fid = fopen(FileN,'w');
fprintf(fid,'%s\n',['# Created on ',datestr(now),' to accompany the ABAQUS Input File: ',NAME.FileName]); 
fprintf(fid,'%s\n',['# (C) Adam Jan Sadowski, 2009']); 
fprintf(fid,'%s\n',' ');

% Writing prerequisites
fprintf(fid,'%s\n','import sys');
fprintf(fid,'%s\n','import os');
fprintf(fid,'%s\n','import visualization');
fprintf(fid,'%s\n','from abaqus import *');
fprintf(fid,'%s\n','from odbAccess import *');
fprintf(fid,'%s\n','from abaqusConstants import *');
fprintf(fid,'%s\n','from odbMaterial import *');
fprintf(fid,'%s\n','from odbSection import *');
fprintf(fid,'%s\n','from caeModules import *');
fprintf(fid,'%s\n','from math import *');
fprintf(fid,'%s\n',' ');

% Writing input requests
fprintf(fid,'%s\n','INC = X  # Write the increment frame (INTEGER) at which the data is desired'); fprintf(fid,'%s\n',' ');

% Opening database
fprintf(fid,'%s\n',['o = session.openOdb(name=''',FileO,''')']);
fprintf(fid,'%s\n','session.viewports[''Viewport: 1''].setValues(displayedObject=o)');
fprintf(fid,'%s\n',' ');

% Cleaning variable names, defensive programming
fprintf(fid,'%s\n','# Initial housekeeping');
for A = 1:length(PATHS.Axial.Generatrix)
    fprintf(fid,'%s\n',['del session.xyDataObjects[''A',num2str(A),'_S22_MEMB'']']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''A',num2str(A),'_S22_BEND'']']);
end
for C = 1:length(PATHS.Circumferential.Generatrix)
    fprintf(fid,'%s\n',['del session.xyDataObjects[''C',num2str(C),'_S22_MEMB'']']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''C',num2str(C),'_S22_BEND'']']);
end
fprintf(fid,'%s\n',' ');

% Writing paths
fprintf(fid,'%s\n','# Axial paths');
for A = 1:length(PATHS.Axial.Generatrix)
    range = [num2str(PATHS.Axial.First(A)),':',num2str(PATHS.Axial.Last(A)),':',num2str(PATHS.Axial.Step(A))];
    fprintf(fid,'%s\n',['session.Path(name=''A',num2str(A),'_Path_',num2str(PATHS.Axial.Generatrix(A)),''', type=NODE_LIST, expression=((''SHELL_INSTANCE'', (''',range,''', )), ))']);
end
fprintf(fid,'%s\n',' '); fprintf(fid,'%s\n','# Circumferential paths');
for C = 1:length(PATHS.Circumferential.Generatrix)
    range = [num2str(PATHS.Circumferential.First(C)),':',num2str(PATHS.Circumferential.Last(C)),':',num2str(PATHS.Circumferential.Step(C))];
    fprintf(fid,'%s\n',['session.Path(name=''C',num2str(C),'_Path_',num2str(PATHS.Circumferential.Generatrix(C)),''', type=NODE_LIST, expression=((''SHELL_INSTANCE'', (''',range,''', )), ))']);
end
fprintf(fid,'%s\n',' ');

% Extracting data requests for S22 (axial) stresses
fprintf(fid,'%s\n','# Extracting data along Axial & Circumferential paths at the required frame');
fprintf(fid,'%s\n','session.viewports[''Viewport: 1''].odbDisplay.setFrame(step=0, frame=INC)');
fprintf(fid,'%s\n','session.viewports[''Viewport: 1''].odbDisplay.setPrimaryVariable(variableLabel=''S'', outputPosition=INTEGRATION_POINT, refinement=(COMPONENT,''S22''))');
fprintf(fid,'%s\n','session.viewports[''Viewport: 1''].odbDisplay.display.setValues(plotState=(CONTOURS_ON_DEF, ))');
for A = 1:length(PATHS.Axial.Generatrix)
    fprintf(fid,'%s\n',['pth = session.paths[''A',num2str(A),'_Path_',num2str(PATHS.Axial.Generatrix(A)),''']']);
    fprintf(fid,'%s\n','session.viewports[''Viewport: 1''].odbDisplay.basicOptions.setValues(sectionResults=USE_BOTTOM)');
    fprintf(fid,'%s\n',['session.XYDataFromPath(name=''A',num2str(A),'_S22_BOT',''', path=pth, includeIntersections=False, shape=UNDEFORMED, labelType=TRUE_DISTANCE_Z)']);
    fprintf(fid,'%s\n',['session.viewports[''Viewport: 1''].odbDisplay.basicOptions.setValues(sectionResults=USE_TOP)']);
    fprintf(fid,'%s\n',['session.XYDataFromPath(name=''A',num2str(A),'_S22_TOP',''', path=pth, includeIntersections=False, shape=UNDEFORMED, labelType=TRUE_DISTANCE_Z)']);
end
for C = 1:length(PATHS.Circumferential.Generatrix)
    fprintf(fid,'%s\n',['pth = session.paths[''C',num2str(C),'_Path_',num2str(PATHS.Circumferential.Generatrix(C)),''']']);
    fprintf(fid,'%s\n','session.viewports[''Viewport: 1''].odbDisplay.basicOptions.setValues(sectionResults=USE_BOTTOM)');
    fprintf(fid,'%s\n',['session.XYDataFromPath(name=''C',num2str(C),'_S22_BOT',''', path=pth, includeIntersections=False, shape=UNDEFORMED, labelType=NORM_DISTANCE)']);
    fprintf(fid,'%s\n',['session.viewports[''Viewport: 1''].odbDisplay.basicOptions.setValues(sectionResults=USE_TOP)']);
    fprintf(fid,'%s\n',['session.XYDataFromPath(name=''C',num2str(C),'_S22_TOP',''', path=pth, includeIntersections=False, shape=UNDEFORMED, labelType=NORM_DISTANCE)']);
end
fprintf(fid,'%s\n',' ');

% Calculating membrane and bending S22 stresses
fprintf(fid,'%s\n','# Calculating membrane and bending stress components');
for A = 1:length(PATHS.Axial.Generatrix)
    fprintf(fid,'%s\n',['xyBot, xyTop = session.xyDataObjects[''A',num2str(A),'_S22_BOT''], session.xyDataObjects[''A',num2str(A),'_S22_TOP'']']);
    fprintf(fid,'%s\n','xyMem, xyBen = 0.5*(xyBot+xyTop), 0.5*(xyBot-xyTop)');
    fprintf(fid,'%s\n',['xyMem.setValues(sourceDescription=''0.5*("A',num2str(A)','_S22_BOT"+"A',num2str(A),'_S22_TOP")'')']);
    fprintf(fid,'%s\n',['xyBen.setValues(sourceDescription=''0.5*("A',num2str(A)','_S22_BOT"-"A',num2str(A),'_S22_TOP")'')']);
    fprintf(fid,'%s\n','tmpNM, tmpNB = xyMem.name, xyBen.name');
    fprintf(fid,'%s\n',['session.xyDataObjects.changeKey(tmpNM, ''A',num2str(A),'_S22_MEMB'')']);
    fprintf(fid,'%s\n',['session.xyDataObjects.changeKey(tmpNB, ''A',num2str(A),'_S22_BEND'')']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''A',num2str(A),'_S22_BOT'']']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''A',num2str(A),'_S22_TOP'']']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''_temp_1'']']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''_temp_3'']']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''_temp_5'']']);
    objB = ['A',num2str(A),'_S22_MEMB']; objT = ['A',num2str(A),'_S22_BEND'];
end
for C = 1:length(PATHS.Circumferential.Generatrix)
    fprintf(fid,'%s\n',['xyBot, xyTop = session.xyDataObjects[''C',num2str(C),'_S22_BOT''], session.xyDataObjects[''C',num2str(C),'_S22_TOP'']']);
    fprintf(fid,'%s\n','xyMem, xyBen = 0.5*(xyBot+xyTop), 0.5*(xyBot-xyTop)');
    fprintf(fid,'%s\n',['xyMem.setValues(sourceDescription=''0.5*("C',num2str(C)','_S22_BOT"+"C',num2str(C),'_S22_TOP")'')']);
    fprintf(fid,'%s\n',['xyBen.setValues(sourceDescription=''0.5*("C',num2str(C)','_S22_BOT"-"C',num2str(C),'_S22_TOP")'')']);
    fprintf(fid,'%s\n','tmpNM, tmpNB = xyMem.name, xyBen.name');
    fprintf(fid,'%s\n',['session.xyDataObjects.changeKey(tmpNM, ''C',num2str(C),'_S22_MEMB'')']);
    fprintf(fid,'%s\n',['session.xyDataObjects.changeKey(tmpNB, ''C',num2str(C),'_S22_BEND'')']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''C',num2str(C),'_S22_BOT'']']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''C',num2str(C),'_S22_TOP'']']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''_temp_1'']']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''_temp_3'']']);
    fprintf(fid,'%s\n',['del session.xyDataObjects[''5'']']);
    objB = ['C',num2str(C),'_S22_MEMB']; objT = ['C',num2str(C),'_S22_BEND'];
end
fprintf(fid,'%s\n',' ');

% Writing to file
fprintf(fid,'%s\n','# Writing to file');
fprintf(fid,'%s\n','xyDATA_A, xyDATA_C = [], []');
for A = 1:length(PATHS.Axial.Generatrix)
    fprintf(fid,'%s\n',['xyDATA_A.append(session.xyDataObjects[''A',num2str(A),'_S22_MEMB''])']);
end
for A = 1:length(PATHS.Axial.Generatrix)
    fprintf(fid,'%s\n',['xyDATA_A.append(session.xyDataObjects[''A',num2str(A),'_S22_BEND''])']);
end
for C = 1:length(PATHS.Circumferential.Generatrix)
    fprintf(fid,'%s\n',['xyDATA_C.append(session.xyDataObjects[''C',num2str(C),'_S22_MEMB''])']);
end
for C = 1:length(PATHS.Circumferential.Generatrix)
    fprintf(fid,'%s\n',['xyDATA_C.append(session.xyDataObjects[''C',num2str(C),'_S22_BEND''])']);
end
fprintf(fid,'%s\n',['session.writeXYReport(fileName=''',FileT_A,''', xyData=xyDATA_A, appendMode=OFF)']);
fprintf(fid,'%s\n',['session.writeXYReport(fileName=''',FileT_C,''', xyData=xyDATA_C, appendMode=OFF)']);
fprintf(fid,'%s\n','print ''2 files written.''');
fclose(fid);
if PATHS.Python.Copy == 1; [a,b,c]=copyfile(FileN,OPCopy.destination,'f'); PATHS.Python.Result = a; end