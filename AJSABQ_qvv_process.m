function [QVV,go] = AJSABQ_qvv_process(QVV,OP)
% Function to process the external .qvv file data
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:24 (previously 28/08/09 - 23:06)

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

go = 1; fide = fopen(OP.cvv_file); kill = 0; J = 0; choice = 0; disp(' '); disp('      Reading data from .qvv file.'); disp(' ');
if fide == -1; kill = 1; go = 0; disp('External .qvv file not found.'); end
while kill == 0
    line = fgetl(fide);
    if line == -1; kill = 1; break; end
    if isempty(line) || length(line) == 1; choice = 0; continue; end
    if length(line) > 1 && choice == 0
        if sum(line(1:2) == '##') == 2; continue; end
        if line(1) == '#' && line(2) ~= '#'; remainder = line(2:length(line));
            switch remainder
                case {'HEADING'} % Choice 1
                    choice = 1; continue
                case {'SOLID_DATA'} % Choice 2
                    choice = 2; continue
                case {'CHANNEL_DATA'} % Choice 3
                    choice = 3; continue
                case {'SILO_DATA'} % Choice 4
                    choice = 4; continue
                case {'PRESSURES'} % Choice 5
                    choice = 5; continue
                case {'END'}
                    kill = 1; break
            end
        end
    end
    if choice > 0
        switch choice
            case {1} % Description
                QVV.Heading = line; continue
            case {2} % Solid data
                if sum(line(1:2) == '##') == 2; continue; end
                num = str2num(line); QVV.SOLID.gam = num(1); QVV.SOLID.mu_i = num(2); QVV.SOLID.mu_w = num(3); 
                QVV.SOLID.Ks_w = num(4); QVV.SOLID.Ks_c = num(5); QVV.SOLID.Kc_w = num(6); QVV.SOLID.Kc_s = num(7); QVV.SOLID.factor = num(8);
            case {3} % Flow channel geometry data
                if sum(line(1:2) == '##') == 2; continue; end
                num = str2num(line); QVV.CHAN.ec = num(1); QVV.CHAN.hJ = num(2); QVV.CHAN.h0 = num(3); QVV.CHAN.power = num(4);
            case {4} % Silo geometry data
                if sum(line(1:2) == '##') == 2; continue; end
                num = str2num(line); H = num(1); R = num(2); Tmax = num(3); QVV.dTT = num(4);
                if R ~= OP.R || H ~= OP.Hc || Tmax ~= OP.T; 
                    go = 0; kill = 1; disp('External .cvv file data mismatch.');
                end % Silo geometry of the .qvv file does not match that of the proposed mesh
            case {5} % Normal pressure data
                if sum(line(1:2) == '##') == 2; continue; end
                J = J + 1; num = str2num(line); QVV.Z(J) = num(1); QVV.T(J) = num(2); QVV.qn(J) = num(3);
        end
    end
end
fclose(fide); disp('      Reading data from .qvv file. Done.'); disp(' ');
for I = 1:size(OP.type,1)
    for J = 1:size(OP.type,2)
        if OP.type(I,J) == 100; OP.factor(I,J) = QVV.SOLID.factor; end
    end
end
