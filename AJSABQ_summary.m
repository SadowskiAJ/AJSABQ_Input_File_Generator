function AJSABQ_summary(TEMP,STEP,RSL,NAME,OP,DRAW,NDS,ELS,GEO,PATHS,SPR,SPRELS,SOLID)
% Function to display a summary of the model in the Matlab workspace 
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:27 (previously 11/02/12 - 16:24)

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

d1 = TEMP.ddd1; dt1 = TEMP.ddt1; d2 = datestr(now); dt2 = clock;
disp(' '); disp(['      Finished writing the input file. ',NAME.FileName]);
W = what; W.path(length(W.path)+1) = '\'; W.path(length(W.path)+1:length(W.path)+length(NAME.FileName)) = NAME.FileName; D = dir(W.path);
if D.bytes/1024 < 1e3; disp(['      Input file size: ',num2str(D.bytes/1024),' kB']);
else disp(['      Input file size: ',num2str(D.bytes/1024e3),' MB']);
end
if STEP.Ana.Restart == 1
    disp(['      Finished writing the restart file. ',NAME.RestartName]);
    WR = what; WR.path(length(WR.path)+1) = '\'; WR.path(length(WR.path)+1:length(WR.path)+length(NAME.RestartName)) = NAME.RestartName; DR = dir(WR.path);
    if DR.bytes/1024 < 1e3; disp(['      Input file size: ',num2str(DR.bytes/1024),' kB']);
    else disp(['      Input file size: ',num2str(DR.bytes/1024e3),' MB']);
    end
end
if OP.Mat_Save == 1; disp('      Input data has been saved to .mat file'); end

if OP.Heltoggle == 1
    if OP.Helorder == 1
        txt = '1st order helical '; 
    elseif OP.Helorder == 2
        txt = '2nd order helical ';
    end
else 
    txt = ''; 
end

disp(' '); 
disp(['      Total nodes: ',num2str(NDS)]); 
disp(['      Total ',txt,'shell elements: ',num2str(ELS)]);
if SPR == 1; disp(['      Total spring elements: ',num2str(SPRELS)]); end

if GEO.feature.Super_ell(1) == 1; disp(' ');
    if GEO.feature.Super_ell(7) == 1; txt = 'Sine on.'; else txt = ''; end
    disp(['      Super-elliptical flattening cylinder feature: inward deflection '...
        ,num2str(GEO.feature.Super_ell(2)),'; tS ',num2str(GEO.feature.Super_ell(3)),'; power p '...
        ,num2str(GEO.feature.Super_ell(4)),'; power q ',num2str(GEO.feature.Super_ell(5)),'; max at height'...
        ,num2str(GEO.feature.Super_ell(6)),'; ',txt]);
end
if GEO.feature.Circ_flat(1) == 1; disp(' ');
    if GEO.feature.Circ_flat(5) == 1; txt = 'Sine on.'; else txt = ''; end
    disp(['      Circular flattening cylinder feature: inward deflection '...
        ,num2str(GEO.feature.Circ_flat(2)),'; angular_spread ',num2str(GEO.feature.Circ_flat(3)),'; max at height'...
        ,num2str(GEO.feature.Circ_flat(4)),'; ',txt]);
end
if GEO.feature.Axi_sin(1) == 1; disp(' ');
    disp(['      Axial sine geometric cylinder feature: ',num2str(GEO.feature.Axi_sin(2)),' inward deflection; '...
        ,num2str(GEO.feature.Axi_sin(3)),' half waves.']);
end
if GEO.feature.Circ_cos(1) == 1; disp(' ');
    if GEO.feature.Circ_cos(4) == 1; txt = 'Sine on.'; else txt = ''; end
    disp(['      Circumferential cosine geometric feature: ',num2str(GEO.feature.Circ_cos(2)),' inward deflection; '...
        ,num2str(GEO.feature.Circ_cos(3)),' waves. ',txt]);
end
if OP.Heltoggle == 1
    if GEO.feature.H_Weld.Toggle == 1 
        disp(' '); disp('      Helical welds included.');       
    end    
end

disp(' ');
for S = 1:STEP.Tot
    s = STEP.Ana.Type(S);
    if s == 1; disp(['      Step ',num2str(S),': Static - General']); end
    if s == 2; disp(['      Step ',num2str(S),': Static - Perturbation']); end
    if s == 3; disp(['      Step ',num2str(S),': Static - Riks']); end
    for T = 1:STEP.Sub(S)
        if OP.factor(S,T) ~= 1; txF = ['; factored by ',num2str(OP.factor(S,T))]; else txF = ''; end
        if OP.ears(S,T) == 0; txE = ['; without ears']; else txE = ''; end
        if OP.EN_Dis(S,T) == 1; txD = ['; with EN 1991-4 concentric discharge factors Ch = ',num2str(SOLID.Ch),' & Cw = ',num2str(SOLID.Cw)]; else txD = ''; end
        switch OP.type(S,T)
            case .1
                disp(['      Concentrated force in X through rigid point (bottom); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case .2
                disp(['      Concentrated force in Y through rigid point (bottom); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);            
            case .3
                disp(['      Concentrated force in Z through rigid point (bottom); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);            
            case .4 
                disp(['      Local moment about X through rigid point (bottom); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case .5
                disp(['      Local moment about Y through rigid point (bottom); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case .6
                disp(['      Local moment about Z through rigid point (bottom); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case .7 
                disp(['      Applied rotation about X through rigid point (bottom); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case .8
                disp(['      Applied rotation about Y through rigid point (bottom); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case .9
                disp(['      Applied rotation about Z through rigid point (bottom); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 1
                disp(['      Concentrated force in X through rigid point (top); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 2
                disp(['      Concentrated force in Y through rigid point (top); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);            
            case 3
                disp(['      Concentrated force in Z through rigid point (top); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);            
            case 4 
                disp(['      Local moment about X through rigid point (top); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 5
                disp(['      Local moment about Y through rigid point (top); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 6
                disp(['      Local moment about Z through rigid point (top); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 7 
                disp(['      Applied rotation about X through rigid point (top); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 8
                disp(['      Applied rotation about Y through rigid point (top); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 9
                disp(['      Applied rotation about Z through rigid point (top); magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 10
                disp(['      Cylinder line load; magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 15
                disp(['      Cylinder internal pressure; magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 152
                disp(['      Cylinder internal friction; magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 16
                disp(['      Hopper internal pressure; magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 162
                disp(['      Hopper internal friction; magnitude = ',num2str(OP.mag(S,T)),txF,'.']);
            case 17
                disp(['      Cylinder internal prescribed distribution; ',OP.Inline.S]);
            case 20
                disp(['      Cylinder Janssen concentric pressures',txF,txD,'.']);
            case 201
                disp(['      Cylinder Janssen concentric pressures + custom inline load',txF,txD,'.']);
            case 25
                disp(['      Cylinder Janssen eccentric pressures; kc = ',num2str(OP.kc(S,T)),'; centered at ',num2str(OP.chan0(S,T)),txF,txE,'.']);
            case 30
                disp(['      Cylinder modified Reimbert concentric pressures',txF,txD,'.']);
            case 35
                disp(['      Cylinder modified Reimbert eccentric pressures; kc = ',num2str(OP.kc(S,T)),'; centered at ',num2str(OP.chan0(S,T)),txF,txE,'.']);
            case 40
                disp(['      Hopper Walker concentric pressures',txF,'.']);
            case 45
                disp(['      Hopper Walker eccentric C1 pressures; eta = ',num2str(OP.hop_eta(S,T)),'; th0 = ',num2str(OP.hop_th0(S,T)),txF,'.']);
            case 46
                disp(['      Hopper Walker eccentric C2 pressures; eta = ',num2str(OP.hop_eta(S,T)),'; th0 = ',num2str(OP.hop_th0(S,T)),txF,'.']);
            case 50
                disp(['      EN 1993-4-1 Appendix C wind pressures (single silo); stagnation = ',num2str(OP.mag(S,T)),txF,'.']);
            case 51
                disp(['      EN 1993-4-1 Appendix C wind pressures (group silo); stagnation = ',num2str(OP.mag(S,T)),txF,'.']);
            case 100
                disp(['      Externally generated .cvv file cylinder pressures',txF,'; file: ',OP.cvv_file]);
        end
    end
end
if OP.CWeld == 1; disp(' '); disp('      Circumferential welds included.'); end
if DRAW.mesh_3D.toggle == 1; disp(' '); disp('      3D Structural mesh drawn.'); end
if ~isempty(DRAW.surf_3D.step); disp('      3D Loading surface drawn.'); end
if PATHS.Toggle == 1; disp(' '); disp('      CAE paths extracted.'); end
if PATHS.Python.Toggle == 1; disp(' '); disp('      Accompanying Python script created.'); end

disp(' '); disp(['      Time elapsed:    ',num2str(etime(dt2,dt1)),' seconds    -    ',num2str(etime(dt2,dt1)/60),' minutes']);

fid = fopen('EXTERNAL_notify.wav');
if fid ~= -1
    [y,Fs] = audioread('EXTERNAL_notify.wav');
    soundsc(y,Fs);
end