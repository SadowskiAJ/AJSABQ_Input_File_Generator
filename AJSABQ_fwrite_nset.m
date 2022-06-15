function [fid] = AJSABQ_fwrite_nset(fid,string,SET,wh,t)
% Function to write the node sets for the .inp file for the 
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:24 (previously 07/06/09 - 13:23)

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

switch wh
    case 'n'
        fprintf(fid,'%s\n',['*NSET, NSET=',string,', INSTANCE=SHELL_INSTANCE']);
    case 'e'
        fprintf(fid,'%s\n',['*ELSET, ELSET=',string,', INSTANCE=SHELL_INSTANCE']);
end        

kill = 0; N = 0; I = 0; node = [];
while kill == 0; N = N + 1; I = I + 1;
    node(N) = SET(I);
    if N == 16
        line = num2str(node); sp = 0;
        for J = 1:length(line)
            if line(J) == ' ' && sp == 0; sp = 1; line(J) = ','; continue; end
            if line(J) == ' ' && sp == 1; continue; end
            if line(J) ~= ' ' && sp == 1; sp = 0; continue; end
            if line(J) ~= ' ' && sp == 0; continue; end
        end
        fprintf(fid,'%s\n',line); node = []; N = 0;
    end
    if I == length(SET) && N ~= 16
        line = num2str(node); sp = 0;
        for J = 1:length(line)
            if line(J) == ' ' && sp == 0; sp = 1; line(J) = ','; continue; end
            if line(J) == ' ' && sp == 1; continue; end
            if line(J) ~= ' ' && sp == 1; sp = 0; continue; end
            if line(J) ~= ' ' && sp == 0; continue; end
        end
        fprintf(fid,'%s\n',line);
        kill = 1;
    end
end
if t == 1; fprintf(fid,'%s\n',['*TRANSFORM, NSET=',string,', TYPE=C']); fprintf(fid,'%1.f, %1.f, %1.f, %1.f, %1.f, %1.f\n',([0,0,0,0,0,1])); end