function [] = AJSABQ_draw_3dsurf(ST_P,ST_F,~,DRAW3DS,S)
% Function to draw a 3D surface through the applied element-by-element loads
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:16 (previously 19/11/09 - 14:59)

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

if (~isempty(ST_P.T) && DRAW3DS == 1) || (~isempty(ST_P.T) && DRAW3DS == 3)
    if ST_P.T(1) == -1 && ST_P.R(1) == -1 && ST_P.Z(1) == -1
    else
        [matT, matZ, matR] = AJSABQ_fvectorise(ST_P.T,ST_P.Z,ST_P.R);
        figure; surf(matT,matR*1e3,matZ,matR*1e3);
        xlabel('Circumferential spread (degrees)'); 
        ylabel('Normal pressure (MPa)'); zlabel('Height above base (mm)');
        title(['Step ',num2str(S),' Normal pressure plot']); view(210,30);
        c = colorbar('Location','SouthOutside'); xlabel(c,'Normal pressures on silo wall, qn [kPa]');
    end
end

if (~isempty(ST_F.T) && DRAW3DS == 2) || (~isempty(ST_F.T) && DRAW3DS == 3)
    if ST_F.T(1) == -1 && ST_F.R(1) == -1 && ST_F.Z(1) == -1
    else
        [matT, matZ, matR] = AJSABQ_fvectorise(ST_F.T,ST_F.Z,ST_F.R);
        figure; surf(matT,matR*1e3,matZ,matR*1e3);
        xlabel('Circumferential spread (degrees)');
        ylabel('Frictional traction (MPa)'); zlabel('Height above base (mm)');
        title(['Step ',num2str(S),' Frictional traction plot']); view(210,30);
        c = colorbar('Location','SouthOutside'); xlabel(c,'Frictional tractions on silo wall, qn [kPa]');
    end
end
