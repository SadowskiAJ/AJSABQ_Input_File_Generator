function [fid] = AJSABQ_write_heading(fid,STEP,NAME,GEO,PP,OP,SOLID,RSL,NODES,ELEMENTSindex,PATHS,SPRING)
% Function to write the Heading of the .inp file
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:40 (previously 06/09/12 - 15:46)

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

disp(' '); disp('Writing heading and introduction.');
fprintf(fid,'%s\n','*HEADING'); 
fprintf(fid,'%s\n',[NAME.Heading,' created on ',datestr(now),' with the ABQAJS input file generator Version 2.0 using ',num2str(PP.CPUs),' processors']);
fprintf(fid,'%s\n','**');
fprintf(fid,'%s\n',['** Shell Max. Radius ',num2str(OP.R),'; Circumferential spread ',num2str(OP.T)]);
if GEO.cyl.Inc == 1
    fprintf(fid,'%s\n',['** Shell Height ',num2str(OP.Hc)]); 
    fprintf(fid,'%s\n',['** ',OP.Class,' shell with H/D ratio ',num2str(OP.Hod)]);
end
if GEO.hop.Inc == 1; fprintf(fid,'%s\n',['** Hopper Height ',num2str(OP.Hc),'; Apex Half-Angle ',num2str(GEO.hop.Beta),' degrees']); end
if GEO.roof.Inc == 1; fprintf(fid,'%s\n',['** Roof Angle ',num2str(GEO.roof.Angle),' degrees']); end

if PATHS.Toggle == 1 && OP.Heltoggle ~= 1
    fprintf(fid,'%s\n','**'); fprintf(fid,'%s\n','** Axial Paths');
    for A = 1:length(PATHS.Axial.Generatrix)
        fprintf(fid,'%s\n',['** Path  ',num2str(PATHS.Axial.First(A)),':',num2str(PATHS.Axial.Last(A)),':',num2str(PATHS.Axial.Step(A)),'  at '...
            ,num2str(PATHS.Axial.Generatrix(A)),' degrees;  ',PATHS.Axial.Description(A,:)]);        
    end
    fprintf(fid,'%s\n','**'); fprintf(fid,'%s\n','** Circumferential Paths');
    for C = 1:length(PATHS.Circumferential.Generatrix)
        if PATHS.Circumferential.Fraction(C,2) == 1 || PATHS.Circumferential.Fraction(C,2) == -1 || PATHS.Circumferential.Fraction(C,2) == 21; txt = 'H(hopper)';
        else txt = 'H(cylinder)'; 
        end
         fprintf(fid,'%s\n',['** Path  ',num2str(PATHS.Circumferential.First(C)),':',num2str(PATHS.Circumferential.Last(C)),':',num2str(PATHS.Circumferential.Step(C)),'  at '...
            ,num2str(PATHS.Circumferential.Generatrix(C)),';  ',num2str(PATHS.Circumferential.Fraction(C,1)),' of ',txt,';  ',PATHS.Circumferential.Description(C,:)]);     
    end
end
 
for S = 1:STEP.Tot
    fprintf(fid,'%s\n','**'); fprintf(fid,'%s\n',['** STEP ',num2str(S)]);
    for I = 1:STEP.Sub(S)
        fprintf(fid,'%s\n',['** LOAD ',num2str(I)]);
        for X = 1:size(STEP.Load.Type,1)
            if STEP.Load.Type(X,1:2) == [S I]
                switch STEP.Load.Type(X,3)
                    case .1
                        fprintf(fid,'%s\n',['** Concentrated force in X through rigid point (bottom); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case .2
                        fprintf(fid,'%s\n',['** Concentrated force in Y through rigid point (bottom); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case .3
                        fprintf(fid,'%s\n',['** Concentrated force in Z through rigid point (bottom); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case .4
                        fprintf(fid,'%s\n',['** Local moment about X through rigid point (bottom); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case .5
                        fprintf(fid,'%s\n',['** Local moment about Y through rigid point (bottom); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case .6
                        fprintf(fid,'%s\n',['** Local moment about Z through rigid point (bottom); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case .7
                        fprintf(fid,'%s\n',['** Applied rotation about X through rigid point (bottom); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case .8
                        fprintf(fid,'%s\n',['** Applied rotation about Y through rigid point (bottom); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case .9
                        fprintf(fid,'%s\n',['** Applied rotation about Z through rigid point (bottom); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 1
                        fprintf(fid,'%s\n',['** Concentrated force in X through rigid point (top); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 2
                        fprintf(fid,'%s\n',['** Concentrated force in Y through rigid point (top); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 3
                        fprintf(fid,'%s\n',['** Concentrated force in Z through rigid point (top); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 4
                        fprintf(fid,'%s\n',['** Local moment about X through rigid point (top); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 5
                        fprintf(fid,'%s\n',['** Local moment about Y through rigid point (top); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 6
                        fprintf(fid,'%s\n',['** Local moment about Z through rigid point (top); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 7
                        fprintf(fid,'%s\n',['** Applied rotation about X through rigid point (top); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 8
                        fprintf(fid,'%s\n',['** Applied rotation about Y through rigid point (top); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 9
                        fprintf(fid,'%s\n',['** Applied rotation about Z through rigid point (top); magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 10
                        fprintf(fid,'%s\n',['** Cylinder top line load; magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 15
                        fprintf(fid,'%s\n',['** Cylinder internal uniform pressure; magnitude = ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 17
                        fprintf(fid,'%s\n',['** Cylinder prescribed internal distribution: ',OP.Inline.S]);
                    case 19
                        if OP.ears(S,I) == 1; txt = ' with ears'; else txt = ''; end
                        fprintf(fid,'%s\n',['** Rotter 1986 eccentric silo pressures centered at theta = ',num2str(OP.chan0(S,I)),...
                            '; kc = ',num2str(OP.kc(S,I)),txt,'; factored by ',num2str(OP.factor(S,I))]);
                    case 20
                        if OP.EN_Dis(S,I) == 1; txEN = '; with EN 1991-4 discharge factors'; else txEN = ''; end
                        fprintf(fid,'%s\n',['** Janssen concentric silo pressures; factored by ',num2str(OP.factor(S,I)),txEN]);
                    case 201
                        if OP.EN_Dis(S,I) == 1; txEN = '; with EN 1991-4 discharge factors'; else txEN = ''; end
                        fprintf(fid,'%s\n',['** Janssen concentric silo pressures + custom inline load; factored by ',num2str(OP.factor(S,I)),txEN]);
                    case 25 
                        if OP.ears(S,I) == 1; txt = ' with ears'; else txt = ''; end
                        fprintf(fid,'%s\n',['** Janssen eccentric silo pressures centered at theta = ',num2str(OP.chan0(S,I)),...
                            '; kc = ',num2str(OP.kc(S,I)),txt,'; factored by ',num2str(OP.factor(S,I))]);
                    case 30
                        if OP.EN_Dis(S,I) == 1; txEN = '; with EN 1991-4 discharge factors'; else txEN = ''; end
                        fprintf(fid,'%s\n',['** Modified Reimbert concentric silo pressures; factored by ',num2str(OP.factor(S,I)),txEN]);
                    case 35
                        if OP.ears(S,I) == 1; txt = ' with ears'; else txt = ''; end
                        fprintf(fid,'%s\n',['** Modified Reimbert eccentric silo pressures centered at theta = ',num2str(OP.chan0(S,I)),...
                            '; kc = ',num2str(OP.kc(S,I)),txt,'; factored by ',num2str(OP.factor(S,I))]);
                    case 40
                        fprintf(fid,'%s\n',['** Walker concentric hopper pressures; factored by ',num2str(OP.factor(S,I))]);
                    case 45
                        fprintf(fid,'%s\n',['** Walker eccentric hopper pressures with smooth depression; theta0 = ',num2str(OP.hop_th0(S,I)),...
                            '; eta = ',num2str(OP.hop_eta(S,I)),'; factored by ',num2str(OP.factor(S,I))]);;
                    case 46
                        fprintf(fid,'%s\n',['** Walker eccentric hopper pressures with hypersmooth depression; theta0 = ',num2str(OP.hop_th0(S,I)),...
                            '; eta = ',num2str(OP.hop_eta(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 50
                        fprintf(fid,'%s\n',['** EN 1993-4-1 Appendix C wind pressures (single silo) with stagnation pressure ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);                        
                    case 51
                        fprintf(fid,'%s\n',['** EN 1993-4-1 Appendix C wind pressures (group silo) with stagnation pressure ',num2str(OP.mag(S,I)),'; factored by ',num2str(OP.factor(S,I))]);
                    case 100
                        fprintf(fid,'%s\n',['** External silo pressure data from .cvv file; factored by ',num2str(OP.factor(S,I))]);
                end
            end
        end
    end
end
fprintf(fid,'%s\n','**');
if GEO.feature.Super_ell(1) == 1
    if GEO.feature.Super_ell(7) == 1; txt = 'Sine on.'; else txt = ''; end
    TXT = ['** Super-elliptical flattening cylinder feature: inward deflection; '...
        ,num2str(GEO.feature.Super_ell(2)),'; theta_S ',num2str(GEO.feature.Super_ell(3)),'; power p '...
        ,num2str(GEO.feature.Super_ell(4)),'; power q ',num2str(GEO.feature.Super_ell(5)),'; max at height '...
        ,num2str(GEO.feature.Super_ell(6)),'; ',txt]; 
    fprintf(fid,'%s\n',TXT); fprintf(fid,'%s\n','**');
end
if GEO.feature.Circ_flat(1) == 1
    if GEO.feature.Circ_flat(5) == 1; txt = 'Sine on.'; else txt = ''; end
    TXT = ['** Circular flattening cylinder feature: inward deflection; '...
        ,num2str(GEO.feature.Circ_flat(2)),'; angular_spread ',num2str(GEO.feature.Circ_flat(3)),'; max at height '...
        ,num2str(GEO.feature.Circ_flat(4)),'; ',txt];
    fprintf(fid,'%s\n',TXT); fprintf(fid,'%s\n','**');
end
if GEO.feature.Axi_sin(1) == 1; disp(' ');
    TXT = ['** Axial sine geometric cylinder feature: ',num2str(GEO.feature.Axi_sin(2)),' inward deflection; '...
        ,num2str(GEO.feature.Axi_sin(3)),' half waves.'];    
    fprintf(fid,'%s\n',TXT); fprintf(fid,'%s\n','**');
end
if GEO.feature.Circ_cos(1) == 1; disp(' ');
    if GEO.feature.Circ_cos(4) == 1; txt = 'Sine on.'; else txt = ''; end
    TXT = ['** Circumferential cosine geometric feature: ',num2str(GEO.feature.Circ_cos(2)),' inward deflection; '...
        ,num2str(GEO.feature.Circ_cos(3)),' waves. ',txt]; 
    fprintf(fid,'%s\n',TXT); fprintf(fid,'%s\n','**');
end
if GEO.feature.Cyl_perturbation(1) == 1; disp(' ');
    TXT = ['** Eigenmode mesh perturbation: n = ',num2str(GEO.feature.Cyl_perturbation(2)),' circumferential full waves, m = '...
        ,num2str(GEO.feature.Cyl_perturbation(3)),' axial half waves, deflection = ',num2str(GEO.feature.Cyl_perturbation(4))];
    fprintf(fid,'%s\n',TXT); fprintf(fid,'%s\n','**');
end
if OP.Sol == 1
    fprintf(fid,'%s\n',['** Solid name; ',SOLID.Name]);
    fprintf(fid,'%s\n',['** Solid wall smoothness type; ',SOLID.Wall]);
    fprintf(fid,'%s\n',['** Solid unit weight = ',num2str(SOLID.Weight),' kN/m3']);
    fprintf(fid,'%s\n',['** Solid lateral pressure ratio (upper) = ',num2str(SOLID.K)]);
    fprintf(fid,'%s\n',['** Solid angle of repose = ',num2str(SOLID.Repose),' degrees']);
    fprintf(fid,'%s\n',['** Solid wall friction coefficient (lower) = ',num2str(SOLID.mewL)]);
    fprintf(fid,'%s\n',['** Solid wall friction coefficient (upper) = ',num2str(SOLID.mewU)]);
    fprintf(fid,'%s\n',['** Internal friction angle (upper) = ',num2str(SOLID.Frang),' degrees']);
    fprintf(fid,'%s\n',['** EN 1991-4 discharge factors: normal, Ch = ',num2str(SOLID.Ch),'; frictional, Cw = ',num2str(SOLID.Cw)]);
    fprintf(fid,'%s\n','**');
end

fprintf(fid,'%s\n',['** RSL.Region = [ ',num2str(RSL.Region),' ]']); 
fprintf(fid,'%s\n',['** RSL.Where = [ ',num2str(RSL.Where),' ]']);
if OP.Heltoggle == 1 
    fprintf(fid,'%s\n',['** RSL.HEL.h1 = [ ',num2str(RSL.HEL.h1),' ]']);
    fprintf(fid,'%s\n',['** RSL.HEL.h2 = [ ',num2str(RSL.HEL.h2),' ]']);     
    fprintf(fid,'%s\n',['** RSL.HEL.Z = [ ',num2str(RSL.HEL.Z),' ]']);
    fprintf(fid,'%s\n',['** RSL.HEL.dp = [ ',num2str(RSL.HEL.dp),' ]']);
    fprintf(fid,'%s\n',['** GEO.helix.element = ',GEO.helix.element]);
    fprintf(fid,'%s\n',['** GEO.helix.thick = ',num2str(GEO.helix.thick)]);
    fprintf(fid,'%s\n',['** GEO.helix.conformal = ',num2str(GEO.helix.conformal)]);
    fprintf(fid,'%s\n',['** GEO.Ttot = ',num2str(GEO.Ttot)]);
    txt = 'helical ';
else
    fprintf(fid,'%s\n',['** RSL.h1 = [ ',num2str(round(RSL.h1)),' ]']);
    fprintf(fid,'%s\n',['** RSL.h2 = [ ',num2str(round(RSL.h2)),' ]']); 
    fprintf(fid,'%s\n',['** RSL.Z = [ ',num2str(RSL.Z),' ]']);
    fprintf(fid,'%s\n',['** RSL.t = [ ',num2str(RSL.t),' ]']);
    fprintf(fid,'%s\n',['** RSL.Type = [ ',RSL.Type,' ]']);
    fprintf(fid,'%s\n',['** RSL.t1 = [ ',num2str(RSL.t1),' ]']);
    fprintf(fid,'%s\n',['** RSL.t2 = [ ',num2str(RSL.t2),' ]']);
    fprintf(fid,'%s\n',['** RSL.T = [ ',num2str(RSL.T),' ]']);
    fprintf(fid,'%s\n',['** RSL.Ttot = [ ',num2str(RSL.Ttot),' ]']);
    txt = '';
end
fprintf(fid,'%s\n','**'); els = [];

fprintf(fid,'%s\n',['** This model contains ',num2str(length(NODES.index)),' nodes']); 
fprintf(fid,'%s\n',['** This model contains ',num2str(length(ELEMENTSindex)),' ',txt,'shell elements']);
if SPRING.toggle == 1
    fprintf(fid,'%s\n',['** This model contains ',num2str(length(SPRING.index)),' spring elements']);
end

if OP.Heltoggle ~= 1
    RSLtdiff = RSL.t2 - RSL.t1; dTmin = min(RSL.t2 - RSL.t1); dTmax = max(RSL.t2 - RSL.t1);
    for E = 1:length(RSL.Type)
        if RSL.Type(E) == 'a'; tip(E,:) = 'S4     '; els(E) = RSL.Z(E)*RSL.Ttot; nZ = RSL.Z(E); nT = RSL.T;
        elseif RSL.Type(E) == 'b'; tip(E,:) = 'S4R    '; els(E) = RSL.Z(E)*RSL.Ttot; nZ = RSL.Z(E); nT = RSL.T;
        elseif RSL.Type(E) == 'c'; tip(E,:) = 'S4R5   '; els(E) = RSL.Z(E)*RSL.Ttot; nZ = RSL.Z(E); nT = RSL.T;
        elseif RSL.Type(E) == 'd'; tip(E,:) = 'S3     '; els(E) = RSL.Z(E)*RSL.Ttot*2; nZ = RSL.Z(E); nT = RSL.T;
        elseif RSL.Type(E) == 'e'; tip(E,:) = 'STRI3  '; els(E) = RSL.Z(E)*RSL.Ttot*2; nZ = RSL.Z(E); nT = RSL.T;
        elseif RSL.Type(E) == 'f'; tip(E,:) = 'STRI65 '; els(E) = RSL.Z(E)*RSL.Ttot/2; nZ = RSL.Z(E)*0.5; nT = RSL.T*0.5;
        elseif RSL.Type(E) == 'g'; tip(E,:) = 'S8R    '; els(E) = RSL.Z(E)*RSL.Ttot/4; nZ = RSL.Z(E)*0.5; nT = RSL.T*0.5;
        elseif RSL.Type(E) == 'h'; tip(E,:) = 'S8R5   '; els(E) = RSL.Z(E)*RSL.Ttot/4; nZ = RSL.Z(E)*0.5; nT = RSL.T*0.5;
        elseif RSL.Type(E) == 'j'; tip(E,:) = 'S9R5   '; els(E) = RSL.Z(E)*RSL.Ttot/4; nZ = RSL.Z(E)*0.5; nT = RSL.T*0.5;
        end
        LZ = (RSL.h2(E) - RSL.h1(E))/nZ; LZ = LZ/sqrt(OP.R*RSL.t(E));
        LT = ((RSL.t2 - RSL.t1)*pi*OP.R/180)./nT; LTmin = min(LT)/sqrt(OP.R*RSL.t(E)); LTmax = max(LT)/sqrt(OP.R*RSL.t(E));
        if length(RSL.t) > 1; ex = [' with thickness t = ',num2str(RSL.t(E)),' mm.']; else ex = '.'; end
        fprintf(fid,'%s\n',['** This input file generates ',num2str(els(E)),' elements of type ',tip(E,:),'between z = ',num2str(RSL.h1(E)),' and ',num2str(RSL.h2(E)),' mm',ex]);
        fprintf(fid,'%s\n',['Axial length: ',num2str(LZ),' sqrt(R*t)']);
        fprintf(fid,'%s\n',['Circ length (min): ',num2str(LTmin),' sqrt(R*t)']);
        fprintf(fid,'%s\n',['Circ length (max): ',num2str(LTmax),' sqrt(R*t)']);
    end
else
    if OP.Heltype == 'a'; tip = 'S4     '; tri = 'S3';
    elseif OP.Heltype == 'b'; tip = 'S4R    '; tri = 'S3';
    elseif OP.Heltype == 'c'; tip = 'S4R5   '; tri = 'S3';
    elseif OP.Heltype == 'g'; tip = 'S8R    '; tri = 'STRI65';
    elseif OP.Heltype == 'h'; tip = 'S8R5   '; tri = 'STRI65';
    elseif OP.Heltype == 'j'; tip = 'S9R5   '; tri = 'STRI65';
    end
    
    fprintf(fid,'%s\n',['** This input file generates ',num2str(OP.Hel_Quad),' QUAD helical elements of type ',tip]);
    fprintf(fid,'%s\n',['** This input file generates ',num2str(OP.Hel_Tri),' TRI helical elements of type ',tri]);
    fprintf(fid,'%s\n',['** QUAD helical element side length #1 (x axis in image plane): ',num2str(OP.D1/sqrt(OP.R*OP.Helt)),' sqrt(R*t)']);
    fprintf(fid,'%s\n',['** QUAD helical element side length #2 (y axis in image plane): ',num2str(OP.D2/sqrt(OP.R*OP.Helt)),' sqrt(R*t)']);
    fprintf(fid,'%s\n',['** QUAD helical element aspect ratio: ',num2str(max(OP.D1/OP.D2,OP.D2/OP.D1))]);
end

if SPRING.toggle == 1
    fprintf(fid,'%s\n',['** This input file generates ',num2str(length(SPRING.index)),' elements of type SPRINGA between z = 0 and ',num2str(OP.Hc),' mm.']);
end

fprintf(fid,'%s\n','**');

if GEO.feature.C_Weld.Toggle == 1
    for W = 1:length(GEO.feature.C_Weld.Z)
        fprintf(fid,'%s\n',['** Circumferential weld of type ',GEO.feature.C_Weld.Type(W),' at z = ',num2str(GEO.feature.C_Weld.Z(W)),' with amplitude ',num2str(GEO.feature.C_Weld.Ampt(W))]);           
    end
    fprintf(fid,'%s\n','**');
end