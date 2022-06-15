function [OP] = AJSABQ_pre_meshing(RSL,OP,PP)
% Function to prepare the meshing allocations
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

i1 = []; i2 = []; e1 = []; e2 = []; h1 = []; h2 = []; Hmin = []; Hmax = []; prt = 0; prh = 0; pre = 0;
for R = 1:length(RSL.Region)
    if RSL.Type(R) == 'g' || RSL.Type(R) == 'h' || RSL.Type(R) == 'j'
        spec = 1; CC = 4; k = 0.25;
    elseif RSL.Type(R) == 'f'; spec = 1; CC = 4; k = 0.5;
    else spec = 0; CC = 2; k = 1;
    end
    if RSL.Type(R) == 'd' || RSL.Type(R) == 'e'; k = 2; end
    RslHtot = RSL.Z(R); div(R,1) = 0; base = RslHtot/PP.CPUs; hdif = RSL.h2(R) - RSL.h1(R); hdiv = hdif/RSL.Z(R); mds = 1;
    while sum(div(R,:)) ~= RslHtot; div(R,:) = []; mds = 1;
        while sum(mds) > 0; mds = [];
            for C = 1:PP.CPUs
                c = ceil(rand(1)*(CC));
                if c == 1; div(R,C) = floor(base);
                elseif c == 2; div(R,C) = ceil(base);
                elseif c == 3; div(R,C) = floor(base)-1;
                elseif c == 4; div(R,C) = ceil(base)+1;
                end
                if CC == 4; mds(C) = mod(div(R,C),2); else mds = []; end
            end
        end
    end
    for D = 1:length(div(R,:))
        if R == 1 && D == 1; i1(R,D) = 1; e1(R,D) = 1; Hmin(R,D) = min(RSL.h1); x = 0; %% Hmin(R,D) = 0
        else i1(R,D) = prt + 1; e1(R,D) = pre + 1; Hmin(R,D) = prh + hdiv; x = 1;
        end
        Hmax(R,D) = Hmin(R,D) + (div(R,D) - x)*hdiv; prh = Hmax(R,D);
        e2(R,D) = e1(R,D) + RSL.Ttot*div(R,D)*k-1; pre = e2(R,D);
        i2(R,D) = i1(R,D) + RSL.Ttot + (div(R,D) - x)*(RSL.Ttot + 1); prt = i2(R,D);
    end
end
OP.i1 = i1; OP.i2 = i2; OP.e1 = e1; OP.e2 = e2; OP.Hmin = Hmin; OP.Hmax = Hmax; OP.div = div;