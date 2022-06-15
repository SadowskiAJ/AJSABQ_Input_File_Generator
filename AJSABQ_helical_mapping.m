function [X,Y, Z] = AJSABQ_helical_mapping(helix,feature,H,R,x,y,NR)
% Helical conformal/non-conformal transformation for the MPI Parallelised helical meshing algorithm
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:24 (previously 06/09/12 - 15:05)

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

adjust = 00;
conformal = helix.conformal;
Pitch = helix.pitch;
alpha = helix.incline*pi/180;
rho = 2*pi*R*sin(alpha);

if conformal == 0
    NT = x; Z = y - Pitch + x*(Pitch/(2*pi));
elseif conformal == 1
    phi = cos(alpha)*(x - y*tan(alpha)); z = y/cos(alpha) + sin(alpha)*(x - y*tan(alpha));
    NT = 2*pi*(phi/(2*pi*R)); Z = z - rho*cos(alpha);
end

% Eigenmode mesh perturbation (if present)
% THE BELOW IS PRETTY HACKY AND AS OF 2022 I DON'T REMEMBER WHAT IT DOES
if feature.Cyl_perturbation(1) == 1
    %NR = NR - feature.Cyl_perturbation(4)*sin(feature.Cyl_perturbation(2)*NT)*sin(0.5*feature.Cyl_perturbation(3)*(Z/H)*(2*pi));

    %     % Temporary alternative
    %     Fourier_Zs = [-0.001059438, -0.745805852, 0.002496493, -0.508437432, 1.133635735, -0.31955094, -0.924383518, -0.222963867, -0.126870896, -0.126689242, 0.016646482, -0.054526238,...
    %         0.049798784, -0.010960989, 0.043399938, 0.008505051, 0.025813068, 0.012291749, 0.010344266, 0.008917276, 0.001405561, 0.004364052, -0.001739362, 0.001343646, -0.001723206,...
    %         0.000119111, -0.000849606, -9.34857E-05, -0.000236359, -3.23756E-05, -2.59669E-05, 1.54725E-06, 4.76677E-07, -3.15937E-07, 1.56806E-07, 8.9212E-08, -3.91361E-07, 1.15409E-07,...
    %         1.4054E-06, -7.48406E-06, -1.39563E-05, -9.75194E-05, -6.02569E-05, -0.000411019, -1.65092E-05, -0.000961072, 0.000502677, -0.001267055, 0.00201881, -0.000120809, 0.004646746,...
    %         0.00412271, 0.007268782, 0.012460928, 0.006865736, 0.023449624, -0.001404409, 0.031005494, -0.023080593, 0.021181368, -0.06260611, -0.039647505, -0.121797201, -0.334795892,...
    %         -0.207696364, 1.006396786, -0.191452544, 0.137497333, -0.36360562, 0.04707476, 0.026574979, 0.312203168, -0.467739967, -0.437036117, -0.339123156, 1.13036779, -0.081616182,...
    %         0.110252563, -0.125352113, 0.005243791, -0.079492413, -0.023677334, -0.038730378, -0.024858204, -0.011562658, -0.016530341, 0.002205812, -0.007663445, 0.006290986, -0.001869571,...
    %         0.005304738, 0.000583571, 0.002910503, 0.00093786, 0.00105195, 0.000542389, 0.000179099, 0.000175514, -3.53462E-05, 2.52552E-05, -2.20767E-05]';
    %
    %     rad = 0; L = length(Fourier_Zs); INTS = [0]; c = 0;
    %     for l = 1:L/2; c = c + 1;
    %         INTS(end+1:end+2) = c*[1 1];
    %     end
    %     if length(INTS) > L; INTS(end) = []; end
    %     if length(INTS) >= 1; coss = [1 2:2:length(INTS)]; else coss = [1]; end
    %     if length(INTS) >= 3; sins = [3:2:length(INTS)]; else sins = []; end
    %
    %     rad = rad + sum(Fourier_Zs(coss).*cos(INTS(coss)'*Z*pi/180));
    %     rad = rad + sum(Fourier_Zs(sins).*sin(INTS(sins)'*Z*pi/180));

    load AJS_spline
    rad = fittedmodel1(Z);
    if Z <= 1e-3 || Z >= H - 1e-3; rad = 0; end

    NR = NR + rad*feature.Cyl_perturbation(4)*sin(feature.Cyl_perturbation(2)*NT);
end

% Angular adjustment
NT = NT + adjust*pi/180;

X = NR*cos(NT); Y = NR*sin(NT);