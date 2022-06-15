function [fid] = AJSABQ_write_material(fid,MATERIAL)
% Function to write the Material of the .inp file
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:40 (previously 11/04/13 - 19:44)

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

disp('Writing material definition.');
fprintf(fid,'%s\n','**');
fprintf(fid,'%s\n',['*MATERIAL, NAME=',MATERIAL.Name]);
if MATERIAL.Type == 1; fprintf(fid,'%s\n','*ELASTIC, TYPE=ISOTROPIC'); fprintf(fid,'%1.f, %2.1f\n',[MATERIAL.E(1);MATERIAL.v(1)]); end
if MATERIAL.Type == 2; fprintf(fid,'%s\n','*ELASTIC, TYPE=LAMINA');
    if length(MATERIAL.G) == 1; MATERIAL.G(2:3) = 0; end
    fprintf(fid,'%1.f, %1.f, %1.1f, %1.f, %1.f, %1.f\n',[MATERIAL.E(1);MATERIAL.E(2);MATERIAL.v(1);MATERIAL.G(1);MATERIAL.G(2);MATERIAL.G(3)]);
end
if isempty(MATERIAL.sy) == 0 || isempty(MATERIAL.ep) == 0
    fprintf(fid,'%s\n','*PLASTIC, HARDENING=ISOTROPIC'); fprintf(fid,'%1.5f, %2.5e\n',[MATERIAL.sy;MATERIAL.ep]);
end
fprintf(fid,'%s\n','**');