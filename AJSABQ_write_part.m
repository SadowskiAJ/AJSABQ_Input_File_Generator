function [fid] = AJSABQ_write_part(fid,RSL,NODES,ELEMENTS,MATERIAL,OP)
% Function to write the Part of the .inp file 
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:40 (previously 28/03/12 - 18:16)

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

Imp = OP.Imp;

disp('Writing part.');
fprintf(fid,'%s\n','*PART, NAME=SHELL_Part');
if OP.Heltoggle == 1
    fprintf(fid,'%s\n','*NODE, SYSTEM=R');
else
    fprintf(fid,'%s\n','*NODE, SYSTEM=C');
end
X = 3; format = ['%1.f, %1.',num2str(X),'f, %1.',num2str(X),'f, %1.',num2str(X),'f\n'];
% Change the 1.Xf to the desired no. of digits after the decimal above
format4 = ['%1.f, %1.f, %1.f, %1.f\n'];
format5 = ['%1.f, %1.f, %1.f, %1.f, %1.f\n'];
format7 = ['%1.f, %1.f, %1.f, %1.f, %1.f, %1.f, %1.f\n'];
format9 = ['%1.f, %1.f, %1.f, %1.f, %1.f, %1.f, %1.f, %1.f, %1.f\n'];
format10 = ['%1.f, %1.f, %1.f, %1.f, %1.f, %1.f, %1.f, %1.f, %1.f, %1.f\n'];
if OP.Heltoggle == 1
    for L = 1:length(NODES.index)
        fprintf(fid,format,[NODES.index(L),NODES.X(L),NODES.Y(L),NODES.Z(L)]);
    end
else
    for L = 1:length(NODES.index)
        fprintf(fid,format,[NODES.index(L),NODES.R(L),NODES.T(L),NODES.Z(L)]);
    end
end

if OP.Heltoggle == 1
    if OP.Heltype == 'a'; tip = 'S4'; frmt = format5; N = 4; end
    if OP.Heltype == 'b'; tip = 'S4R'; frmt = format5; N = 4; end
    if OP.Heltype == 'c'; tip = 'S4R5'; frmt = format5; N = 4; end
    if OP.Heltype == 'g'; tip = 'S8R'; frmt = format9; N = 8; end
    if OP.Heltype == 'h'; tip = 'S8R5'; frmt = format9; N = 8; end
    if OP.Heltype == 'j'; tip = 'S9R5'; frmt = format10; N = 9; end
    
    fprintf(fid,'%s\n','**'); fprintf(fid,'%s\n',['*ELEMENT, TYPE=',tip]);
    for L = 1:OP.Hel_Quad
        fprintf(fid,frmt,[ELEMENTS.index(L),ELEMENTS.nodes(L,(1:N))]);     
    end
    
    if OP.Helorder == 1; tip = 'S3'; frmt = format4; N = 3; end
    if OP.Helorder == 2; tip = 'STRI65'; frmt = format7; N = 6; end
    fprintf(fid,'%s\n','**'); fprintf(fid,'%s\n',['*ELEMENT, TYPE=',tip]);      
    for L = (OP.Hel_Quad+1):(OP.Hel_Quad+OP.Hel_Tri)
        fprintf(fid,frmt,[ELEMENTS.index(L),ELEMENTS.nodes(L,(1:N))]);           
    end    
    
else
    for E = 1:length(RSL.Region)
        if char(RSL.Type(E)) == 'a'; tip = 'S4'; frmt = format5; N = 4; end % 4-node general purpose, double curved, finite membrane strains (suitable for large-strain analysis), thick and thin shells
        if char(RSL.Type(E)) == 'b'; tip = 'S4R'; frmt = format5; N = 4; end % 4-node reduced integration, hourglass control, finite membrane strains
        if char(RSL.Type(E)) == 'c'; tip = 'S4R5'; frmt = format5; N = 4; end % S4R with small membrane strains (5 degrees of freedom per node), thin shells only
        if char(RSL.Type(E)) == 'd'; tip = 'S3'; frmt = format4; N = 3; end % 3-node triangular general purpose, double curved, finite membrane strains, thick and thin shells
        if char(RSL.Type(E)) == 'e'; tip = 'STRI3'; frmt = format4; N = 3; end % 3-node triangular facet thin shell element
        if char(RSL.Type(E)) == 'f'; tip = 'STRI65'; frmt = format7; N = 6; end % 6-node triangular (5 degrees of freedom per node), thin shells only
        if char(RSL.Type(E)) == 'g'; tip = 'S8R'; frmt = format9; N = 8; end % 8-node general purpose, double curved, thick shell, small membrane strains, reduced integration
        if char(RSL.Type(E)) == 'h'; tip = 'S8R5'; frmt = format9; N = 8; end % S8R with 5 degrees of freedom per node, thin shells only
        if char(RSL.Type(E)) == 'j'; tip = 'S9R5'; frmt = format10; N = 9; end % 9-node with 5 degrees of freedom per node, thin shells only
        fprintf(fid,'%s\n','**'); fprintf(fid,'%s\n',['*ELEMENT, TYPE=',tip]);
        for L = 1:length(ELEMENTS.region)
            if ELEMENTS.region(L) == RSL.Region(E)
                fprintf(fid,frmt,[ELEMENTS.index(L),ELEMENTS.nodes(L,(1:N))]);
            end
        end
    end
end

fprintf(fid,'%s\n','**');
if OP.Heltoggle == 1
    fprintf(fid,'%s\n','*NSET, NSET=NTOTAL, GENERATE');
    fprintf(fid,'%1.f, %1.f, %1.f\n',[min(NODES.index),max(NODES.index),1]);
    fprintf(fid,'%s\n','*ELSET, ELSET=ETOTAL, GENERATE');
    fprintf(fid,'%1.f, %1.f, %1.f\n',[min(ELEMENTS.index),max(ELEMENTS.index),1]);
else
    fprintf(fid,'%s\n','*NSET, NSET=NTOTAL, GENERATE');
    fprintf(fid,'%1.f, %1.f, %1.f\n',[min(NODES.index),max(NODES.index),1]);
    fprintf(fid,'%s\n','*ELSET, ELSET=ETOTAL, GENERATE');
    fprintf(fid,'%1.f, %1.f, %1.f\n',[min(ELEMENTS.index),max(ELEMENTS.index),1]);
end

if OP.Heltoggle == 1
    fprintf(fid,'%s\n','**');
    fprintf(fid,'%s\n',['*SHELL SECTION, ELSET=ETOTAL, MATERIAL=',MATERIAL.Name,', SECTION INTEGRATION=SIMPSON, POISSON=ELASTIC']); 
    fprintf(fid,'%s\n',[num2str(OP.Helt),', 99']);
else
    tsofar = [];
    for T = 1:length(RSL.t); tsofar(T) = RSL.t(T); clear('ex'); fprintf(fid,'%s\n','**');
        if isempty(find(tsofar == RSL.t(T), 1)) == 0; ex = ['_',num2str(length(find(tsofar==RSL.t(T))))]; else ex = (''); end
        if T == 1; ebot = 1; etop = length(find(ELEMENTS.region == T));
        else
            ebot = 1 + find(ELEMENTS.region == (T-1), 1, 'last' ); etop = find(ELEMENTS.region == (T-1), 1, 'last' ) + length(find(ELEMENTS.region == T));
        end
        t = num2str(RSL.t(T)); a = find(t == '.'); t(a) = 'p';
        fprintf(fid,'%s\n',['*ELSET, ELSET=EL_',t,ex,', GENERATE']); fprintf(fid,'%1.f, %1.f, %1.f\n',[ebot,etop,1]);
        fprintf(fid,'%s\n',['*SHELL SECTION, ELSET=EL_',t,ex,', MATERIAL=',MATERIAL.Name,', SECTION INTEGRATION=SIMPSON, POISSON=ELASTIC']); fprintf(fid,'%s\n',[num2str(RSL.t(T)),', 99']);
    end
end
fprintf(fid,'%s\n','*END PART'); fprintf(fid,'%s\n','**');

if Imp.toggle == 1
    %fprintf(fid,'%s\n',['*IMPERFECTION, FILE=',Imp.file,', STEP=1, INC=X']);
    fprintf(fid,'%s\n',['*IMPERFECTION, FILE=',Imp.file,', STEP=',num2str(OP.Steps)]);
    fprintf(fid,'%s\n',['1, ',num2str(Imp.ampt)]); fprintf(fid,'%s\n','**');
end