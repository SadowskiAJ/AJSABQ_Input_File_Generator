function [fid,NSET] = AJSABQ_write_assembly(fid,GEO,OP,RSL,STEP,NSET,ELSET,NODES,ELEMENTS,SPRING)
% Function to write the Assembly of the .inp file
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:27 (previously 24/01/13 - 19:44)

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

disp('Writing assembly and node/element sets.');
fprintf(fid,'%s\n','**');
fprintf(fid,'%s\n','*ASSEMBLY, NAME=SHELL_ASSEMBLY');
fprintf(fid,'%s\n','*INSTANCE, NAME=SHELL_INSTANCE, PART=SHELL_PART');
fprintf(fid,'%s\n','*END INSTANCE');
fprintf(fid,'%s\n','**');

if ~isempty(find(STEP.BC.Top_RigidRing(:,2:7) == 1, 1))
    fprintf(fid,'%s\n','*NODE');
    fprintf(fid,'%1.f, %1.f, %1.f, %1.f\n',[1, 0, 0, OP.H]);
    fprintf(fid,'%s\n','*NSET, NSET=RIGID_NODE_TOP');
    fprintf(fid,'%s\n',num2str(1));
%     if OP.Heltoggle == 1;
%         fprintf(fid,'%s\n',['*TRANSFORM, NSET=RIGID_NODE_TOP, TYPE=C']); 
%         fprintf(fid,'%s\n',['0, 0, 0, 0, 0, 1']);
%     end 
    fprintf(fid,'%s\n','**');
end
if ~isempty(find(STEP.BC.Trans_RigidRing(:,2:7) == 1, 1))
    fprintf(fid,'%s\n','*NODE');
    fprintf(fid,'%1.f, %1.f, %1.f, %1.f\n',[2, 0, 0, 0]);
    fprintf(fid,'%s\n','*NSET, NSET=RIGID_NODE_BOTTOM');
    fprintf(fid,'%s\n',num2str(2));
   % if OP.Heltoggle == 1;
%         fprintf(fid,'%s\n',['*TRANSFORM, NSET=RIGID_NODE_BOTTOM, TYPE=C']); 
%         fprintf(fid,'%s\n',['0, 0, 0, 0, 0, 1']);
   % end 
end
if SPRING.toggle == 1 && SPRING.Discharge ~= 1
    fprintf(fid,'%s\n','**'); fprintf(fid,'%s\n',['*ELEMENT, TYPE=SPRINGA, ELSET=SPRINGS']);
    for S = 1:length(SPRING.index)
        fprintf(fid,'%s\n',[num2str(SPRING.index(S)),', SHELL_INSTANCE.',num2str(SPRING.nodes(S,1)),', SHELL_INSTANCE.',num2str(SPRING.nodes(S,2))]);
    end
    fprintf(fid,'%s\n','*SPRING, ELSET=SPRINGS'); fprintf(fid,'%s\n',''); fprintf(fid,'%s\n',num2str(SPRING.constant));
    fprintf(fid,'%s\n','**');
end
if SPRING.toggle == 1 && SPRING.Discharge == 1; SET = 1; cont = 0;
    for S = 1:length(SPRING.set)
        if SPRING.set(S) == SET && cont == 0; cont = 1;
            fprintf(fid,'%s\n','**'); fprintf(fid,'%s\n',['*ELEMENT, TYPE=SPRINGA, ELSET=SPRINGS_',num2str(SPRING.setlist.no(SET))]);
            fprintf(fid,'%s\n',[num2str(SPRING.index(S)),', SHELL_INSTANCE.',num2str(SPRING.nodes(S,1)),', SHELL_INSTANCE.',num2str(SPRING.nodes(S,2))]);
        elseif SPRING.set(S) == SET && cont == 1
            fprintf(fid,'%s\n',[num2str(SPRING.index(S)),', SHELL_INSTANCE.',num2str(SPRING.nodes(S,1)),', SHELL_INSTANCE.',num2str(SPRING.nodes(S,2))]);
        elseif SPRING.set(S) ~= SET && cont == 1; cont = 2;
            fprintf(fid,'%s\n',['*SPRING, ELSET=SPRINGS_',num2str(SPRING.setlist.no(SET))]); fprintf(fid,'%s\n',''); fprintf(fid,'%s\n',num2str(SPRING.setlist.k(SET)));
        end
        if SPRING.set(S) == (SET + 1) && cont == 2; cont = 1; SET = SET + 1;
            fprintf(fid,'%s\n','**'); fprintf(fid,'%s\n',['*ELEMENT, TYPE=SPRINGA, ELSET=SPRINGS_',num2str(SPRING.setlist.no(SET))]);
            fprintf(fid,'%s\n',[num2str(SPRING.index(S)),', SHELL_INSTANCE.',num2str(SPRING.nodes(S,1)),', SHELL_INSTANCE.',num2str(SPRING.nodes(S,2))]);
        end
    end
    fprintf(fid,'%s\n',['*SPRING, ELSET=SPRINGS_',num2str(SPRING.setlist.no(SET))]); fprintf(fid,'%s\n',''); fprintf(fid,'%s\n',num2str(SPRING.setlist.k(SET)));
end

if ~isempty(find(STEP.Load.Type(:,3) == 10, 1)) % Cylinder top edge line load
    if OP.Heltoggle ~= 1
        fprintf(fid,'%s\n','*ELSET, ELSET=ETOP, INSTANCE=SHELL_INSTANCE, GENERATE');
        if RSL.Where(length(RSL.Where)) == 3; Eltop = RSL.Type(length(RSL.Type)-1); Ecum = sum(RSL.Els(1:length(RSL.Els)-1));
        elseif RSL.Where(length(RSL.Where)) == 2; Eltop = RSL.Type(length(RSL.Type)); Ecum = sum(RSL.Els(1:length(RSL.Els)));
        end
        if Eltop == 'a' || Eltop == 'b' || Eltop == 'c'; Rsltop = RSL.Ttot; step = 1; adj = 1; txt = 'ETOP, E3'; end
        if Eltop == 'd' || Eltop == 'e'; Rsltop = RSL.Ttot*2; step = 2; adj = 2; txt = 'ETOP, E1'; end
        if Eltop == 'f'; Rsltop = RSL.Ttot; step = 2; adj = 2; txt = 'ETOP, E1'; end
        if Eltop == 'g' || Eltop == 'h' || Eltop == 'j'; Rsltop = RSL.Ttot/2; step = 1; adj = 1; txt = 'ETOP, E3'; end
        fprintf(fid,'%1.f, %1.f, %1.f\n',[Ecum - Rsltop + adj,Ecum,step]);
        fprintf(fid,'%s\n','*SURFACE, TYPE=ELEMENT, NAME=TOP_EDGE');
        fprintf(fid,'%s\n',txt); fprintf(fid,'%s\n','**'); % For 4 and 8-node elements, E3 is the element surface corresponding to the top edge, for 3 and 6-node elements, it is E1.
    elseif OP.Heltoggle == 1
        [fid] = AJSABQ_fwrite_nset(fid,'ETOP',ELSET.Top,'e',0); 
        fprintf(fid,'%s\n','*SURFACE, TYPE=ELEMENT, NAME=TOP_EDGE');
        fprintf(fid,'%s\n','ETOP, E3');
    end
end

if ~isempty(find(STEP.Load.Type(:,3) == 15, 1)) || ~isempty(find(STEP.Load.Type(:,3) == 152, 1)); Ecum = 0; E1 = 0; % Cylinder internal pressure/friction
    for L = 1:length(RSL.Where); Ecum = Ecum + RSL.Els(L);
        if RSL.Where(L) == 2 && E1 == 0; E1 = Ecum - RSL.Els(L) + 1; end
        if RSL.Where(L) == 3; E2 = Ecum - RSL.Els(L); end
        if RSL.Where(L) == 2 && L == length(RSL.Where); E2 = Ecum; end
    end
    fprintf(fid,'%s\n','*ELSET, ELSET=CYL_ETOTAL, INSTANCE=SHELL_INSTANCE, GENERATE');
    fprintf(fid,'%1.f, %1.f, %1.f\n',[E1,E2,1]);
    fprintf(fid,'%s\n','*SURFACE, TYPE=ELEMENT, NAME=INNER_SURFACE_CYL, INTERNAL');
    fprintf(fid,'%s\n','CYL_ETOTAL, SNEG'); fprintf(fid,'%s\n','**');
end

if ~isempty(find(STEP.Load.Type(:,3) == 16, 1)) || ~isempty(find(STEP.Load.Type(:,3) == 162, 1)); Ecum = 0; E1 = 1; % Hopper internal pressure/friction
    for L = 1:length(RSL.Where); Ecum = Ecum + RSL.Els(L);
        if RSL.Where(L) == 2 || L == length(RSL.Where); E2 = Ecum - RSL.Els(L); break; end
    end
    fprintf(fid,'%s\n','*ELSET, ELSET=HOP_ETOTAL, INSTANCE=SHELL_INSTANCE, GENERATE');
    fprintf(fid,'%1.f, %1.f, %1.f\n',[E1,E2,1]);
    fprintf(fid,'%s\n','*SURFACE, TYPE=ELEMENT, NAME=INNER_SURFACE_HOP, INTERNAL');
    fprintf(fid,'%s\n','HOP_ETOTAL, SNEG'); fprintf(fid,'%s\n','**');
end

if OP.Heltoggle ~= 1
    if ~isempty(find(STEP.BC.T0(:,2) == 1, 1)) || ~isempty(find(STEP.BC.T0(:,2) == -1, 1)) ||...
            ~isempty(find(STEP.BC.T0(:,2) == 2, 1)) || ~isempty(find(STEP.BC.T0(:,2) == -2, 1))
        [fid] = AJSABQ_fwrite_nset(fid,'LEFT',NSET.Left,'n',0); fprintf(fid,'%s\n','**');
    end
    if ~isempty(find(STEP.BC.Ttot(:,2) == 1, 1)) || ~isempty(find(STEP.BC.Ttot(:,2) == -1, 1)) ||...
            ~isempty(find(STEP.BC.Ttot(:,2) == 2, 1)) || ~isempty(find(STEP.BC.Ttot(:,2) == -2, 1))
        [fid] = AJSABQ_fwrite_nset(fid,'RIGHT',NSET.Right,'n',0); fprintf(fid,'%s\n','**');
    end
    if GEO.cyl.Inc == 1
        if ~isempty(find(STEP.BC.Top(:,2:7) == 1, 1)) || ~isempty(find(STEP.BC.Top_RigidRing(:,2:7) == 1, 1))
            NSET.Top = NSET.Top(length(NSET.Top)-length(NSET.Bottom)+1:length(NSET.Top));
            if ~isempty(STEP.BC.Top_RigidRing); t = 0; else t = 1; end
            [fid] = AJSABQ_fwrite_nset(fid,'TOP',NSET.Top,'n',t); fprintf(fid,'%s\n','**');
        end
    end
end

if ~isempty(find(STEP.BC.Trans(:,2:7) == 1, 1)) || ~isempty(find(STEP.BC.Trans_RigidRing(:,2:7) == 1, 1))
    if ~isempty(STEP.BC.Top_RigidRing); t = 0; else t = 1; end
    [fid] = AJSABQ_fwrite_nset(fid,'TRANSITION',NSET.Bottom,'n',t); 
    %if OP.Heltoggle == 1;
        fprintf(fid,'%s\n',['*TRANSFORM, NSET=TRANSITION, TYPE=C']); 
        fprintf(fid,'%s\n',['0, 0, 0, 0, 0, 1']);
        fprintf(fid,'%s\n','**');
        [fid] = AJSABQ_fwrite_nset(fid,'TOP',NSET.Top,'n',t); 
        fprintf(fid,'%s\n',['*TRANSFORM, NSET=TOP, TYPE=C']); 
        fprintf(fid,'%s\n',['0, 0, 0, 0, 0, 1']);
        fprintf(fid,'%s\n','**'); 
    %end
end

if GEO.hop.Inc == 1
    if ~isempty(find(STEP.BC.Pit(:,2:7) == 1, 1)); [fid] = AJSABQ_fwrite_nset(fid,'PIT',NSET.Czelusc,'n',0); end
    if ~isempty(find(STEP.BC.Pit(:,8) == 1, 1)); fprintf(fid,'%s\n','**');
        for N = 1:length(NSET.Czelusc) % transformation into local coordinate sys for every one of the bottom nodes of the hopper
            fprintf(fid,'%s\n',['*NSET, NSET=PIT_',num2str(N),', INSTANCE=SHELL_INSTANCE']);
            fprintf(fid,'%s\n',[num2str(N)]);
            fprintf(fid,'%s\n',['*TRANSFORM, NSET=PIT_',num2str(N),', TYPE=C']);
            Ra = NODES.R(N); Ta = NODES.T(N); Za = NODES.Z(N);
            Rb = NODES.R(N+RSL.Ttot); Tb = NODES.T(N+RSL.Ttot); Zb = NODES.Z(N+RSL.Ttot);
            xa = Ra*cos(Ta*pi/180); ya = Ra*sin(Ta*pi/180); za = Za; xb = Rb*cos(Tb*pi/180); yb = Rb*sin(Tb*pi/180); zb = Zb;
            fprintf(fid,'%1.2f, %1.2f, %1.2f, %1.2f, %1.2f, %1.2f\n',([xa,ya,za,xb,yb,zb])); fprintf(fid,'%s\n','**');
        end
    else
        fprintf(fid,'%s\n',['*TRANSFORM, NSET=PIT, TYPE=C']); fprintf(fid,'%1.f, %1.f, %1.f, %1.f, %1.f, %1.f\n',([0,0,0,0,0,1])); fprintf(fid,'%s\n','**');
    end
end

if SPRING.toggle == 1
    [fid] = AJSABQ_fwrite_nset(fid,'AXIS',NSET.Axis,'n',0); fprintf(fid,'%s\n','**');
end
if OP.Riks == 1
    NSET.Riks = NSET.Top(1);
    if length(NSET.Riks) > 1; NSET.Riks(2:length(NSET.Riks)) = []; end
    if OP.Heltoggle ~= 1
        Nopp = find(round(NODES.Z) == OP.H & NODES.T == OP.T);
        NSET.Path.Centre = [num2str(NSET.Riks),':1:-',num2str(RSL.Ttot)];
        NSET.Path.Oppopsite = [num2str(Nopp),':',num2str(RSL.Ttot + 1),':-',num2str(RSL.Ttot)];
    end    
    fprintf(fid,'%s\n',['*NSET, NSET=RIKS, INSTANCE=SHELL_INSTANCE']); fprintf(fid,'%s\n',num2str(NSET.Riks)); fprintf(fid,'%s\n','**');
end
if OP.Heltoggle == 1
    [fid] = AJSABQ_fwrite_nset(fid,'N_LEFTS',OP.P.Lefts.ID,'n',0); fprintf(fid,'%s\n','**');
    [fid] = AJSABQ_fwrite_nset(fid,'N_RIGHTS',OP.P.Rights.ID,'n',0); fprintf(fid,'%s\n','**');
    [fid] = AJSABQ_fwrite_nset(fid,'N_INTS',OP.P.Ints.ID,'n',0); fprintf(fid,'%s\n','**');
    [fid] = AJSABQ_fwrite_nset(fid,'N_TOPS',OP.P.Tops.ID,'n',0); fprintf(fid,'%s\n','**');
    if GEO.helix.order == 2
        [fid] = AJSABQ_fwrite_nset(fid,'N_LEFTS2',OP.P.Lefts2.ID,'n',0); fprintf(fid,'%s\n','**');
        [fid] = AJSABQ_fwrite_nset(fid,'N_RIGHTS2',OP.P.Rights2.ID,'n',0); fprintf(fid,'%s\n','**');
        [fid] = AJSABQ_fwrite_nset(fid,'N_INTS2',OP.P.Ints2.ID,'n',0); fprintf(fid,'%s\n','**');
        [fid] = AJSABQ_fwrite_nset(fid,'N_TOPS2',OP.P.Tops2.ID,'n',0); fprintf(fid,'%s\n','**');
    end
end

if ~isempty(find(STEP.BC.Top_RigidRing(:,2:7) == 1, 1))
    fprintf(fid,'%s\n','*SURFACE, TYPE=NODE, NAME=TOP_SURFACE');
    fprintf(fid,'%s\n','TOP');
    fprintf(fid,'%s\n','*COUPLING, CONSTRAINT NAME=RIGID_TOP_RING, REF NODE=RIGID_NODE_TOP, SURFACE=TOP_SURFACE');
    fprintf(fid,'%s\n','*KINEMATIC'); fprintf(fid,'%s\n','**');
end
if ~isempty(find(STEP.BC.Trans_RigidRing(:,2:7) == 1, 1))
    fprintf(fid,'%s\n','*SURFACE, TYPE=NODE, NAME=BOTTOM_SURFACE');
    fprintf(fid,'%s\n','TRANSITION');
    fprintf(fid,'%s\n','*COUPLING, CONSTRAINT NAME=RIGID_BOTTOM_RING, REF NODE=RIGID_NODE_BOTTOM, SURFACE=BOTTOM_SURFACE');
    fprintf(fid,'%s\n','*KINEMATIC'); fprintf(fid,'%s\n','**');
end

fprintf(fid,'%s\n','*END ASSEMBLY'); fprintf(fid,'%s\n','**');