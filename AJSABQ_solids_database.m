function [Weight,K,Repose,h0,mewU,mewL,Frang,Ch,Cw] = AJSABQ_solids_database(Name,Type,Rtot)
% Function database holding material properties from EN 1994-1 Table E.1
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:24 (previously 17/12/09 - 10:34)

% Weight - unit weight - kN/m3 (upper characteristic value)
% K - lateral pressure ratio - dimensionless (upper characteristic value)
% Repose - angle of repose - degrees
% h0 - value of z at the highest solid-wall contact - modified Reimbert distributions only
% mewU - wall friction coefficient - dimensionless (upper characteristic value)
% mewL - wall friction coefficient - dimensionless (lower characteristic value, required for flow channel - EN 1991-4:2006 Section 5.2.4.3.1)
% Frang - internal angle of friction - degrees (upper characteristic value)

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

Ch = 1.15; % discharge factor for horizontal wall pressures - EN 1991-4
Cw = 1.1; % discharge factor for wall frictional tractions  - EN 1991-4

switch Name
    case 'Wheat'
        Weight = 9; % unit weight - kN/m3 (upper characteristic value)
        K = 0.54; aK = 1.11; % lateral pressure ratio - dimensionless (upper characteristic value)
        Repose = 34; % angle of repose - degrees
        switch Type % wall friction coefficient - dimensionless
            case 'D1'
                mew = 0.24;
            case 'D2'
                mew = 0.38;
            case 'D3'
                mew = 0.57;             
        end
        aM = 1.16;
        Frang = 30; aF = 1.12; % internal angle of friction - degrees
    case 'Cement'
        Weight = 16; % unit weight - kN/m3 (upper characteristic value)
        K = 0.54; aK = 1.2; % lateral pressure ratio - dimensionless (upper characteristic value)
        Repose = 36; % angle of repose - degrees
        switch Type % wall friction coefficient - dimensionless
            case 'D1'
                mew = 0.41;
            case 'D2'
                mew = 0.46;
            case 'D3'
                mew = 0.51;             
        end
        aM = 1.07;
        Frang = 30; aF = 1.22; % internal angle of friction - degrees  
end

K = K*aK; 
h0 = Rtot*tan(Repose*pi/180)/3;
mewU = mew*aM; 
mewL = mew/aM;
Frang = Frang*aF;