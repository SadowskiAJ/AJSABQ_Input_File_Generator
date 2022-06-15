function [] = AJSABQ_draw_3dmesh(ELEMENTS,NODES,SPRING,OP,DRAW_OP,Ttot)
% Function to draw the structural mesh of the shell
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 14:19 (previously 28/02/12 - 17:06)

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

% NOTE - Element types; a - S4; b - S4R; c - S4R5; d - S3; e - STRI3; f - STRI65; g - S8R; h - S8R5; j - S9R5
figure; hold on;
ELEMENTS.type = char(ELEMENTS.type); SN = 0;
if SPRING.toggle == 1
   for S = 2:length(SPRING.nodes(:,1)); if SPRING.nodes(S,1) ~= SPRING.nodes(S-1,1); SN = SN + 1; end; end
   SN = SN + 1; Sold = 0;
end
for E = 1:length(ELEMENTS.index)
    x = []; y = []; z = [];
    if ELEMENTS.type(E) == 'a' || ELEMENTS.type(E) == 'b' || ELEMENTS.type(E) == 'c'; N = 4; end
    if ELEMENTS.type(E) == 'd' || ELEMENTS.type(E) == 'e'; N = 3; end
    if ELEMENTS.type(E) == 'f'; N = 6; end
    if ELEMENTS.type(E) == 'g' || ELEMENTS.type(E) == 'h'; N = 8; end
    if ELEMENTS.type(E) == 'j'; N = 9; end
    for D = 1:N; EL = ELEMENTS.nodes(E,D);
        if OP.Roof > 0 && ELEMENTS.nodes(E,D) > length(NODES.R)-SN; EL = length(NODES.R)-SN; end
        if OP.Heltoggle == 1
            x(D) = NODES.X(find(NODES.index == EL));
            y(D) = NODES.Y(find(NODES.index == EL));
            z(D) = NODES.Z(find(NODES.index == EL));
        elseif OP.Heltoggle ~= 1
            x(D) = NODES.R(EL)*cos(NODES.T(EL)*pi/180);
            y(D) = NODES.R(EL)*sin(NODES.T(EL)*pi/180);
            z(D) = NODES.Z(EL);
        end
    end
    if N == 3 || N == 6
        fill3([x(1) x(2) x(3) x(1)],[y(1) y(2) y(3) y(1)],[z(1) z(2) z(3) z(1)],'g');
        if DRAW_OP.nodes == 1
            for n = 1:3
                text(x(n),y(n),z(n),num2str(ELEMENTS.nodes(E,n)),'FontWeight','bold','FontSize',8,'Color',[1 0 0]);
            end
        end
    end
    if N == 4 || N == 8 || N == 9
        fill3([x(1) x(2) x(3) x(4) x(1)],[y(1) y(2) y(3) y(4) y(1)],[z(1) z(2) z(3) z(4) z(1)],'r');
        if DRAW_OP.nodes == 1
            for n = 1:4
                text(x(n),y(n),z(n),num2str(ELEMENTS.nodes(E,n)),'FontWeight','bold','FontSize',8,'Color',[1 0 0]);
            end
        end
    end
    if DRAW_OP.elements == 1
        xM = mean([min(x) max(x)]); yM = mean([min(y) max(y)]); zM = mean([min(z) max(z)]);
        text(xM,yM,zM,num2str(E),'FontWeight','bold','FontSize',8);
    end
end
if SPRING.toggle == 1
    for S = 1:length(SPRING.index)
        if OP.Roof > 0; dif = Ttot; else dif = 0; end
        i1 = SPRING.nodes(S,1) - dif; i2 = SPRING.nodes(S,2);  
        x1 = NODES.R(i1)*cos(NODES.T(i1)*pi/180); x2 = NODES.R(i2)*cos(NODES.T(i2)*pi/180);
        y1 = NODES.R(i1)*sin(NODES.T(i1)*pi/180); y2 = NODES.R(i2)*sin(NODES.T(i2)*pi/180);
        z1 = NODES.Z(i1); z2 = NODES.Z(i2);
        line([x1 x2],[y1 y2],[z1 z2],'Color',[137 157 164]/255,'LineStyle','--'); Snew = i1;
        if DRAW_OP.nodes == 1 && Snew ~= Sold 
            text(x1,y1,z1,num2str(i1),'FontWeight','bold','FontSize',8,'Color',[1 0 0]); Sold = Snew;
        end
        if DRAW_OP.elements == 1
            xM = mean([x1 x2]); yM = mean([y1 y2]); zM = mean([z1 z2]);
            text(xM,yM,zM,['S',num2str(S)],'FontWeight','bold','FontSize',8,'Color','b');
        end
    end
end
axis([-OP.R OP.R -OP.R OP.R 0 OP.Hc]); view(70,30); xlabel('X coordinate'); ylabel('Y coordinate'); zlabel('Z coordinate'); grid on;
