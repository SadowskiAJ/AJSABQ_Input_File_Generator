% ABAQUS QVV File Generator Version 1.0
% Numerical solution of ODE based on slice theory to generate the
% granular solid unsymmetrical pressure distributions in a silo
% IMPORTANT - the theoretical background to these equations has since been
% revisited and improved, but the equations coded here reflect those of Dr
% Sadowski's PhD thesis
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 15:24 (previously 30/04/10 - 00:43)

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

clear all, clc

H = 14; % m - total height of silo (no external hopper)
R = 3.4; % m - silo radius
dZ = H/1e3; % m - vertical step size
dT = 1; % degrees - circumferential step size
Tmax = 180; % degrees - max. circumferential extent of 3D surface plot
factor = 0.5;  % factor the pressures
draw = [6]; scribe = 0; File = 'AJS_Thesis_B_silo_ECCec100'; 
Scribe = 'Thesis concentric model, B silo; H=14, R=3.4, ec=(almost100%), n=1.5, r0=0.25';
%File = 'AJS_Thesis_B_silo_pressures';
% 1 - line drawing of vertical pressures
% 11 - line drawing of normal pressures
% 2 - 2D contour plot of vertical pressures
% 3 - area, perimeter and angular components
% 4 - 3D cylindrical surface plot of normal pressures
% 5 - 3D straight axes surface plot of normal pressures
% 6 - 3D geometry plot of silo and flow channel
request = [0 45 90 135 180]; % Starting from the left, 0 <= ANG <= 180
sammanfattning = 1; % request summary output of data to Matlab window
ID = 0;
if ID == 0
ec = 0.0*R; % m - eccentricity of flow channel from the radius
n = 1.5; % power of the distribution, must be >= 1; n = 1 is for linear (conical flow channel)
r0 = 0.25; % m - radius of channel at outlet (z = H)
% Initial study
elseif ID == 1; ec = 0; n = 10; r0 = 0.25; % CON_P
elseif ID == 101; ec = 0; n = 5; r0 = 0.25; % CON_P shallower
elseif ID == 2; ec = 0; n = 1.2; r0 = 0.25; % CON_M
elseif ID == 3; ec = 1; n = 1.2; r0 = 0.25; % ECC_M
elseif ID == 4; ec = 2.5; n = 2; r0 = 0.25; % ECC_P
% Power study
elseif ID == 50;ec = 2.72; n = 5.0; r0 = 0.25; % ECCn5p0
elseif ID == 5; ec = 2.72; n = 4.5; r0 = 0.25; % ECCn4p5
elseif ID == 6; ec = 2.72; n = 4.0; r0 = 0.25; % ECCn4p0
elseif ID == 7; ec = 2.72; n = 3.5; r0 = 0.25; % ECCn3p5
elseif ID == 8; ec = 2.72; n = 3.0; r0 = 0.25; % ECCn3p0
elseif ID == 9; ec = 2.72; n = 2.5; r0 = 0.25; % ECCn2p5
elseif ID == 10; ec = 2.72; n= 2.0; r0 = 0.25; % ECCn2p0
elseif ID == 11; ec = 2.72; n= 1.5; r0 = 0.25; % ECCn1p5
elseif ID == 12; ec = 2.72; n=1.05; r0 = 0.25; % ECCn1p0
% Eccentricity study
elseif ID == 13; ec = 1.00*R-0.25; n = 1.5; r0 = 0.25; % ECCec100
elseif ID == 14; ec = 0.75*R; n = 1.5; r0 = 0.25; % ECCec075
elseif ID == 15; ec = 0.50*R; n = 1.5; r0 = 0.25; % ECCec050
elseif ID == 16; ec = 0.25*R; n = 1.5; r0 = 0.25; % ECCec025
elseif ID == 17; ec = 0.00*R; n = 1.5; r0 = 0.25; % ECCec000
end

MAT = 1;
if MAT == 1 % Wheat
gam = 9000; % N/m3 - unit weight of granular solid
mu_i = 0.664398; % internal friction coefficient
mu_w = 0.4408; % wall friction coefficient
Ch = 1; % normal discharge coefficient
elseif MAT == 2 % Cement
   gam = 16000; mu_i = 0.742665541; mu_w = 0.4922; Ch = 1;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Hereafter please do not modify the code unless you know exactly what you are doing %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a_i = atan(mu_i); a_w = atan(mu_w); % friction coefficients in terms of angles
At = pi*R*R; % m2 - cross-sectional area of silo
bCR = atan(sqrt(tan(a_i)^2+1)-tan(a_i)); % radians - critical internal channel slope angle (betaCR)
h0 = r0/(n*tan(bCR)); % m - virtual origin of parabolic channel, below the base of the silo 
m = r0/(h0^(1/n)); % expansion 'coefficient' of the parabolic channel equation
z12 = H + h0 - ((R-ec)/m)^n; % m - Region 1-2 boundary depth
z23 = H + h0 - ((R+ec)/m)^n; % m - Region 2-3 boundary depth
b12 = atan(m*(H + h0 - z12)^(1/n)/(n*(H + h0 - z12)))*180/pi; % degrees - Channel slope angle at region 1-2 boundary
b23 = atan(m*(H + h0 - z23)^(1/n)/(n*(H + h0 - z23)))*180/pi; % degrees - Channel slope angle at region 1-2 boundary

Kas_i = (1+sin(a_i)^2-2*sqrt(sin(a_i)^2-mu_i^2*cos(a_i)^2))/(4*mu_i^2+cos(a_i)^2); % AS3774 sliding & plastic failure value - internal
Kas_w = (1+sin(a_i)^2-2*sqrt(sin(a_i)^2-mu_w^2*cos(a_i)^2))/(4*mu_w^2+cos(a_i)^2); % AS3774 sliding & plastic failure value - vs. wall
Ks_w = Kas_w; % Static material vs. wall 
Ks_c = Kas_i; % Static material vs. flowing channel 
Kc_w = Kas_w; % Flowing channel vs. wall 
Kc_s = Kas_i; % Flowing channel vs. static material

qC = []; qS = []; qJ = []; I = 0; f = 1; T = 0; qW = []; 
Z = []; rSL = []; rSR = []; rC = []; Zz = []; qrC = []; qrSL = []; qrSR = []; tx = struct('tx',[],'z',[]); 
srf = struct('X',[],'Y',[],'Z',[],'T',[],'P',[]); srf_I = 0; sI = 0; lC = 0;
ch = struct('X',[],'Y',[],'Z',[]); st = struct('X',[],'Y',[],'Z',[]); lc = struct('x',[],'y',[],'z',[]);
EN = struct('pC',[],'pS',[]); ROT = struct('pC',[],'pS',[]);
Fz = []; AS = []; AC = []; UWC = []; UWS = []; USC = []; A_T = []; A_D = [];
% Uwc - contact perimeter between silo wall and flow channel
% Uws - contact perimeter between silo wall and stationary solid
% Usc - contact perimeter between stationary solid and flow channel
% Ac - cross-sectional area of flow channel
% As - cross-sectional area of stationary solid
% ThetaC - flow channel wall contact angle w.r.t silo centre
% DeltaC - flow channel wall contact angle w.r.t channel centre

if z23 <= 0 && z12 <= 0; range_Z = [0:(dZ/1e0):(H-dZ) (H-dZ):(dZ/1e2):H]; 
elseif z23 <= 0 && z12 > 0; range_Z = [0:(dZ/1e0):(z12-dZ) (z12-dZ):(dZ/1e2):(z12+dZ) (z12+dZ):(dZ/1e0):(H-dZ) (H-dZ):(dZ/1e2):H];
elseif z23 > 0 && z12 > 0 && z12 > z23; range_Z = [0:(dZ/1e0):(z23-dZ) (z23-dZ):(dZ/1e2):(z23+dZ)...
        (z23+dZ):(dZ/1e0):(z12-dZ) (z12-dZ):(dZ/1e2):(z12+dZ) (z12+dZ):(dZ/1e0):(H-dZ) (H-dZ):(dZ/1e2):H];
elseif z12 == z23 && z12 > 0; range_Z = [0:(dZ/1e0):(z12-dZ) (z12-dZ):(dZ/1e2):(z12+dZ) (z12+dZ):(dZ/1e0):(H-dZ) (H-dZ):(dZ/1e2):H];
end
range_Z = [0:dZ:H];
range_T = [0:dT:Tmax];
if scribe == 1; File(length(File)+1:length(File)+4) = '.qvv'; if isempty(find(draw == 4, 1)); draw(length(draw)+1) = 4; end
    % Note: because the input file generator assumes mm and N everywhere, the transformation will be done here
    fid = fopen(File,'w'); 
    fprintf(fid,'%s\n',['## Created on ',datestr(now)]); fprintf(fid,'%s\n',['## (C) Adam Jan Sadowski, 2010']); fprintf(fid,'%s\n',' ');
    fprintf(fid,'%s\n','#HEADING'); fprintf(fid,'%s\n',Scribe); fprintf(fid,'%s\n',' '); 
    fprintf(fid,'%s\n','#SOLID_DATA'); fprintf(fid,'%s\n','## gamma (N/mm3), mu_i (-), mu_w (-), Ks_wall (-), Ks_chan (-), Kc_wall (-), Kc_static (-), factor (-)'); 
    fprintf(fid,'%s\n',[num2str(gam/1e9),', ',num2str(mu_i),', ',num2str(mu_w),', ',num2str(Ks_w),', ',num2str(Ks_c),', ',num2str(Kc_w),', ',num2str(Kc_s),', ',num2str(factor)]); 
    fprintf(fid,'%s\n',' ');
    fprintf(fid,'%s\n','#CHANNEL_DATA'); fprintf(fid,'%s\n','## ec (mm), z12 (mm), z23 (mm), power (-), chan_rad (mm), h0 (mm), betaCR (degrees)'); 
    fprintf(fid,'%s\n',[num2str(ec*1e3),', ',num2str(z12*1e3),', ',num2str(z23*1e3),', ',num2str(n),', ',num2str(r0*1e3),', ',num2str(h0*1e3),', ',num2str(bCR*180/pi)]); 
    fprintf(fid,'%s\n',' ');
    fprintf(fid,'%s\n','#SILO_DATA'); fprintf(fid,'%s\n','## H (mm), R (mm), Tmax (degrees), dT (degrees)');
    fprintf(fid,'%s\n',[num2str(H*1e3),', ',num2str(R*1e3),', ',num2str(Tmax),', ',num2str(dT)]); fprintf(fid,'%s\n',' ');
    fprintf(fid,'%s\n','#PRESSURES'); fprintf(fid,'%s\n','## Wall normal pressures in a cylindrical CS: Z (mm), T (degrees), qn (MPa - N/mm2)'); 
end

if sammanfattning == 1; clc; % Output of data summary to the Matlab command window
    disp(' '); disp(' Silo'); disp([' Height: ',num2str(H),' m']); disp([' Radius: ',num2str(R),' m']); disp([' Spread: ',num2str(Tmax),' degrees']);
    disp(' '); disp(' Channel'); disp([' Eccentricity: ',num2str(ec),' m']); disp([' Outlet radius: ',num2str(r0),' m']); disp([' Power: ',num2str(n)]);
    disp([' Beta CR: ',num2str(bCR*180/pi),' degrees']); disp([' Virtual offset: ',num2str(h0),' m']); disp([' Expansion coefficient: ',num2str(m)]);
    if z12 >= 0; disp([' Region 1-2 depth: ',num2str(z12),' m']); end
    if z23 >= 0; disp([' Region 2-3 depth: ',num2str(z23),' m']); end
    if z12 >= 0; disp([' Region 1-2 slope angle: ',num2str(b12),' degrees']); end
    if z23 >= 0; disp([' Region 2-3 slope angle: ',num2str(b23),' degrees']); end
    disp(' '); disp(' Solid'); disp([' Unit weight: ',num2str(gam/1e3),' kN/m3']); disp([' Internal friction angle: ',num2str(a_i*180/pi),' degrees']);
    disp([' Wall friction coefficient: ',num2str(mu_w)]); disp([' K AS3774 vs wall: ',num2str(Kas_w)]); disp([' K AS3774 vs itself: ',num2str(Kas_i)]);
    disp(' ');
end

for I = 1:length(range_Z)-1
    z = range_Z(I); Z(I) = z; dz = range_Z(I+1)-range_Z(I);
    if I == 1; qC(I) = 0; qS(I) = 0; qJ(I) = 0; else qC(I) = qC(I-1) + dqC; qS(I) = qS(I-1) + dqS; qJ(I) = qJ(I-1) + dqJ; dqC = 0; dqS = 0; dqJ = 0; end
    rc1 = m*(H + h0 - z)^(1/n); rc2 = m*(H + h0 - (z + dz))^(1/n); rcm = 0.5*(rc1 + rc2);
    beta1 = atan(m*(H + h0 - z)^(1/n)/(n*(H + h0 - z)));
    beta2 = atan(m*(H + h0 - (z + dz))^(1/n)/(n*(H + h0 - (z + dz)))); beta = 0.5*(beta1 + beta2);
    c = cos(beta); s = sin(beta); t = tan(beta); cot = 1/t;
    Fz(I) = (t*t - Kas_i)/(t*t + 2*mu_i*t - 1); % internal flow channel wall pressure ratio    
    dqJ = (gam - qJ(I)/(At/(mu_w*Kc_w*2*pi*R)))*dz;
    
    if z >= 0 && z < z23
        % This is Region 3 (optional if channel forms at the surface) with mass flow and axisymmetric Janssen pressures; qS = 0 here always.
        Uwc = 2*pi*R; Uws = 0; Usc = 0; Ac = pi*R*R; zo = At/(mu_w*Kc_w*Uwc);
        dqC = (gam - qC(I)/zo)*dz; dqS = 0;
        rC(I,1) = 0; rC(I,2) = 2*R; rSL(I,1) = 0; rSL(I,2) = 0; rSR(I,1) = 2*R; rSR(I,2) = 2*R;
        if isempty(find(tx.tx == 3, 1)) == 1; T = T + 1; tx.tx(T) = 3; tx.z(T) = H - 0.5*z23; end
        UWC(I) = Uwc; UWS(I) = Uws; USC(I) = Usc; AC(I) = Ac; AS(I) = 0; Tc1 = 0; A_T(I) = 180; A_D(I) = 180;
        
    elseif z >= z23
        if R + ec - rc1 < 0; rC(I,1) = 0; rSL(I,1) = 0; rSL(I,2) = 0; else rC(I,1) = R + ec - rc1; rSL(I,1) = 0; rSL(I,2) = rC(I,1); end
        if R + ec + rc1 > 2*R; rC(I,2) = 2*R; rSR(I,1) = 2*R; rSR(I,2) = 2*R; else rC(I,2) = R + ec + rc1; rSR(I,1) = rC(I,2); rSR(I,2) = 2*R; end
        if z >= z23 && z < z12
            % This is Region 2 (optional if channel is eccentric, and one part forms agains the silo wall).
            Tc1 = acos((R*R+ec*ec-rc1*rc1)/(2*R*ec));
            Tc2 = acos((R*R+ec*ec-rc2*rc2)/(2*R*ec));
            if Tc1 <= pi/2 && asin(R*sin(Tc2)/rc2) < asin(R*sin(Tc1)/rc1); D1 = asin(R*sin(Tc1)/rc1); else D1 = pi - asin(R*sin(Tc1)/rc1); end
            Uwc = 2*Tc1*R; Uws = 2*(pi-Tc1)*R; Usc = 2*(pi-D1)*rc1;
            Ac1 = (pi-D1)*rc1*rc1+Tc1*R*R-R*ec*sin(Tc1);
            if isreal(Tc2) == 1
                if Tc2 <= pi/2 && asin(R*sin(Tc2)/rc2) < asin(R*sin(Tc1)/rc1); D2 = asin(R*sin(Tc2)/rc2); else D2 = pi - asin(R*sin(Tc2)/rc2); end
                Ac2 = (pi-D2)*rc2*rc2+Tc2*R*R-R*ec*sin(Tc2);
                lC = lC + 1; lc.x(lC) = R*cos(Tc2); lc.y(lC) = R*sin(Tc2); lc.z(lC) = H-z;
            else
                Tc2 = 0; D2 = 0; Ac2 = pi*rc2*rc2;
            end
            if z23 < 0; tp = 0; else tp = z23; end
            if isempty(find(tx.tx == 2, 1)) == 1; T = T + 1; tx.tx(T) = 2; tx.z(T) = H-(0.5*(z12-tp)+tp); end
                        
        elseif z >= z12 && z <= H
            % This is Region 1 (always present) and represents a fully internal flow channel
            Tc1 = 0; Tc2 = 0; D1 = 0; D2 = 0;
            Uwc = 0; Uws = 2*pi*R; Usc = 2*pi*rc1;
            Ac1 = pi*rc1*rc1; Ac2 = pi*rc2*rc2;
            if isempty(find(tx.tx == 1, 1)) == 1; T = T + 1; tx.tx(T) = 1; tx.z(T) = H-(0.5*(H-z12)+z12); end
        end
        
        As1 = At - Ac1; As2 = At - Ac2; 
        dqC = qC(I)*(Ac1-Ac2)/Ac2 + gam*(Ac1+Ac2)*dz/(2*Ac2) - qC(I)*(mu_w*Kc_w*Uwc+Fz(I)*(t+mu_i)*Usc)*dz/Ac2; 
        if f == 1; qS(I) = qC(I)*Fz(I)*(t+mu_i)/(t+mu_w*Ks_w); qS0 = qS(I); qC0 = qC(I); qJ0 = qJ(I); f = 0; end
        dqS = qS(I)*(As1-As2)/As2 + gam*(As1+As2)*dz/(2*As2) - qS(I)*mu_w*Ks_w*Uws*dz/As2 + qC(I)*Fz(I)*(t+mu_i)*Usc*dz/As2;
        UWC(I) = Uwc; UWS(I) = Uws; USC(I) = Usc; AC(I) = Ac1; AS(I) = As1; A_T(I) = 0.5*(Tc1+Tc2)*180/pi; A_D(I) = 0.5*(D1+D2)*180/pi;
    end
  
    for A = 1:length(request)
        if z >= 0 && z < z23; qW(I,A) = qC(I)*Kc_w/1e3; end % Region 3
        if z >= z23 && z < z12 % Region 2
            if request(A) <= Tc1*180/pi || request(A) >= (360 - Tc1*180/pi); qW(I,A) = qC(I)*Kc_w/1e3;
            elseif request(A) >= Tc1*180/pi && request(A) <= (360 - Tc1*180/pi); qW(I,A) = qS(I)*Ks_w/1e3;
            end
        end
        if z >= z12 && z <= H; qW(I,A) = qS(I)*Ks_w/1e3; end % Region 1
    end
    
    if isempty(find(draw == 4, 1)) == 0 || isempty(find(draw == 5, 1)) == 0 || isempty(find(draw == 6, 1)) == 0
        %ch = struct('X',[],'Y',[],'Z',[]);
        srf_I = srf_I + 1; srf_J = 0; sI = sI + 1; sJ = 0;
        for t = 1:length(range_T); srf_J = srf_J + 1;
            srf.X(srf_I,srf_J) = R*cos(range_T(t)*pi/180); srf.Y(srf_I,srf_J) = R*sin(range_T(t)*pi/180); 
            srf.T(srf_I,srf_J) = range_T(t); srf.Z(srf_I,srf_J) = H-z; 
            if z >= 0 && z < z23 % Region 3
                srf.P(srf_I,srf_J) = qC(I)*Kc_w; 
                sJ = sJ + 1; ch.X(sI,sJ) = R*cos(range_T(t)*pi/180); ch.Y(sI,sJ) = R*sin(range_T(t)*pi/180); ch.Z(sI,sJ) = H-z;
                st.X(sI,sJ) = R*cos(range_T(t)*pi/180); st.Y(sI,sJ) = R*sin(range_T(t)*pi/180); st.Z(sI,sJ) = H-z;
            end 
            if z >= z23 && z < z12 % Region 2
                if range_T(t) <= Tc1*180/pi || range_T(t) >= (360 - Tc1*180/pi) % Inside the perimeter contact
                    srf.P(srf_I,srf_J) = qC(I)*Kc_w;
                    sJ = sJ + 1; ch.X(sI,sJ) = R*cos(range_T(t)*pi/180); ch.Y(sI,sJ) = R*sin(range_T(t)*pi/180); ch.Z(sI,sJ) = H-z;
                elseif range_T(t) >= Tc1*180/pi && range_T(t) <= (360 - Tc1*180/pi) % Outside the perimeter contact
                    srf.P(srf_I,srf_J) = qS(I)*Ks_w;
                    sJ = sJ + 1; ch.X(sI,sJ) = ec + rcm*cos(range_T(t)*pi/180); ch.Y(sI,sJ) = rcm*sin(range_T(t)*pi/180); ch.Z(sI,sJ) = H-z;
                    if sqrt(ch.X(sI,sJ)^2 + ch.Y(sI,sJ)^2) >= R; ch.X(sI,sJ) = ch.X(sI,sJ-1); ch.Y(sI,sJ) = ch.Y(sI,sJ-1); end
                end
                st.X(sI,sJ) = R*cos(range_T(t)*pi/180); st.Y(sI,sJ) = R*sin(range_T(t)*pi/180); st.Z(sI,sJ) = H-z;                
            end
            if z >= z12 && z <= H % Region 1
                srf.P(srf_I,srf_J) = qS(I)*Ks_w; 
                sJ = sJ + 1; ch.X(sI,sJ) = ec + rcm*cos(range_T(t)*pi/180); ch.Y(sI,sJ) = rcm*sin(range_T(t)*pi/180); ch.Z(sI,sJ) = H-z;
                st.X(sI,sJ) = R*cos(range_T(t)*pi/180); st.Y(sI,sJ) = R*sin(range_T(t)*pi/180); st.Z(sI,sJ) = H-z;
            end 
            if scribe == 1
                fprintf(fid,'%s\n',[num2str((H-z)*1e3),', ',num2str(range_T(t)),', ',num2str(srf.P(srf_I,srf_J)*factor/1e6)]);
            end
        end
    end

    qrC(I,1) = qC(I)/1000; qrC(I,2) = qC(I)/1000; 
    qrSL(I,1) = qS(I)/1000; qrSL(I,2) = qS(I)/1000; 
    qrSR(I,1) = qS(I)/1000; qrSR(I,2) = qS(I)/1000; 
    Zz(I,1) = H-Z(I); Zz(I,2) = H-Z(I);
end
Z = H - Z; qC = qC/1000; qS = qS/1000; qJ = qJ/1000;  mx = max(max(max(qC*Ch,qS*Ch),max(qC*Ch,qJ*Ch)));

if sammanfattning == 1
    I12 = find(Z == H - z12); I23 = find(Z == H - z23); IE = length(Z); disp(' Vertical Pressures');
    ICm = find(qC == max(qC*Ch)); ISm = find(qS == max(qS*Ch)); IJm = find(qJ == max(qJ*Ch));
    if z12 >= 0 % Boundary 1/2 pressures
        disp([' Boundary 1/2 pressures @ y/H = ',num2str(1-z12/H)]); disp([' Stationary: ',num2str(qS(I12)*Ch),' kPa']); disp([' Channel: ',num2str(qC(I12)*Ch),' kPa']); disp([' Janssen: ',num2str(qJ(I12)*Ch),' kPa']); disp(' ');
    end
    if z23 >= 0 % Boundary 2/3 pressures
        disp([' Boundary 2/3 pressures @ y/H = ',num2str(1-z23/H)]); disp([' Stationary: ',num2str(qS(I23)*Ch),' kPa']); disp([' Channel: ',num2str(qC(I23)*Ch),' kPa']); disp([' Janssen: ',num2str(qJ(I23)*Ch),' kPa']); disp(' ');
    end
    % Base pressures
    disp(' Base pressures'); disp([' Stationary: ',num2str(qS(IE)*Ch),' kPa']); disp([' Channel: ',num2str(qC(IE)*Ch),' kPa']); disp([' Janssen: ',num2str(qJ(IE)*Ch),' kPa']); disp(' ');
    % Max. pressures
    disp(' Max. pressures'); disp([' Stationary: ',num2str(qS(ISm)*Ch),' kPa @ y/H = ',num2str(Z(ISm)/H)]);
    disp([' Channel: ',num2str(qC(ICm)*Ch),' kPa @ y/H = ',num2str(Z(ICm)/H)]); disp([' Janssen: ',num2str(qJ(IJm)*Ch),' kPa @ y/H = ',num2str(Z(IJm)/H)]); disp(' ');
end
disp('Finished integrating.');
colours = [[1 0.6941 0.3922]; [0 0.498 0]; [1 0 0]; [0 0 1]; [1 0 1]; [1 .502 .502]; [.502 0 .502]; [.502 .502 0]]; % Orange, Green, Red, Blue, Purple, Pinkish, Dark Purple, Greenish
if isempty(find(draw == 1, 1)) == 0 % Line plots - VERTICAL pressure components 
    disp('Plotting lines.'); figure1 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
    axes1 = axes('FontSize',14,'Parent',figure1);
    plot(qC*Ch,Z,'Color',[1 0.6941 0.3922],'LineWidth',3); hold on; plot(qS*Ch,Z,'Color',[0 0.498 0],'LineWidth',3); plot(qJ*Ch,Z,'r','LineWidth',3);
    xlabel('Solid vertical pressures, pv (kPa)','FontSize',14); ylabel('Vertical coordinate, y (m)','FontSize',14); grid on; axis(axes1,[0 mx 0 H]);
    plot([0,mx],[H-z23,H-z23],'k','LineWidth',3,'LineStyle','--'); plot([0,mx],[H-z12,H-z12],'k','LineWidth',3,'LineStyle',':');
    txt1 = ['Flowing Channel, dz = ',num2str(dZ)]; txt2 = ['Stationary Solid, dz = ',num2str(dZ)]; txt3 = ['Janssen Solid, dz = ',num2str(dZ)]; 
    TXT = {txt1,txt2,txt3};
    if z23 >= 0; txt4 = ['Region 3/2, z23 = ',num2str(H-z23)]; TXT(length(TXT)+1) = {txt4}; end
    if z12 >= 0; txt5 = ['Region 2/1, z12 = ',num2str(H-z12)]; TXT(length(TXT)+1) = {txt5}; end
    legend1 = legend(axes1,TXT,'FontSize',12);
end
if isempty(find(draw == 11, 1)) == 0 % Line plots - NORMAL pressure components
    disp('Plotting requested lines.'); figure11 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
    axes1 = axes('FontSize',14,'Parent',figure11); mx = max(max(qW*Ch));
    for A = 1:length(request)
        n = length(num2str(request(A))); 
        if n == 1; no = '  '; elseif n == 2; no = ' '; elseif n == 3; no = ''; end
        plot(qW(:,A)*Ch,Z,'Color',colours(A,:),'LineWidth',3); hold on; 
        txt(A,1:11) = ['Angle = ',no,num2str(request(A))];
    end
    xlabel('Solid normal pressures, pn (kPa)','FontSize',14); ylabel('Vertical coordinate, y (m)','FontSize',14); grid on; axis(axes1,[0 mx 0 H]);
    plot([0,mx],[H-z23,H-z23],'k','LineWidth',3,'LineStyle','--'); plot([0,mx],[H-z12,H-z12],'k','LineWidth',3,'LineStyle',':');
    %txt1 = ['Region 3/2, z23 = ',num2str(H-z23)]; txt2 = ['Region 2/1, z12 = ',num2str(H-z12)];
    legend1 = legend(axes1,{txt},'FontSize',12);
end
if isempty(find(draw == 2, 1)) == 0; Font = 16; % Filled contour pressure plot - 2D
    levels = 100; disp('Plotting 2D contours.');
    figure2 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
    axes2 = axes('FontSize',Font,'Layer','top','Parent',figure2);
    contourf(rC,Zz,qrC*Ch,levels,'LineStyle','none'); hold on; contourf(rSL,Zz,qrSL*Ch,levels,'LineStyle','none'); contourf(rSR,Zz,qrSR*Ch,levels,'LineStyle','none');
    %xlabel('r - radial coordinate'); 
    ylabel('y - vertical coordinate (m)'); set(axes2,'YGrid','on'); offset_R = 2*R; offset_Z = 0; axis([-offset_R 2*R+offset_R -offset_Z H+offset_Z]);
    c = colorbar('Location','SouthOutside','Box','on','FontSize',Font); xlabel(c,'Vertical pressure in solid, pv [kPa]');
    txt = ['H/D = ',num2str(H/(2*R)),'; dz = H/',num2str(H/dZ),'; ec = ',num2str(ec),'; m = ',num2str(m)]; TXT = {}; clear('txt3','txt4');
    TXT = {'Flow Channel','Stationary Solid (Left)','Stationary Solid (Right)'};
    if z23 >= 0; txt3 = ['Region 3/2, y23 = ',num2str(H-z23)]; TXT(length(TXT)+1) = {txt3}; end
    if z12 >= 0; txt4 = ['Region 2/1, y12 = ',num2str(H-z12)]; TXT(length(TXT)+1) = {txt4}; end
    plot([-offset_R,2*R+offset_R],[H-z23,H-z23],'k','LineWidth',3,'LineStyle','--'); 
    plot([-offset_R,2*R+offset_R],[H-z12,H-z12],'k','LineWidth',3,'LineStyle',':');
    title(txt); legend(TXT);
    for t = 1:T; text(-1.5*R,tx.z(t),['Region ',num2str(tx.tx(t))]); end
end
if isempty(find(draw == 3, 1)) == 0 % Line plots - geometric components
    disp('Plotting geometric components.');
    figure; plot(AC,Z,'b','LineWidth',3); hold on; plot(AS,Z,'r','LineWidth',3); grid on;
    xlabel('Area components (m2)'); ylabel('y - vertical coordinate (m)'); legend('Area of Channel - Ac','Area of Stationary Solid - As');
    figure; plot(UWC,Z,'b','LineWidth',3); hold on; plot(UWS,Z,'r','LineWidth',3); plot(USC,Z,'Color',[0 0.498 0],'LineWidth',3); grid on;
    xlabel('Perimeter components (m)'); ylabel('z - vertical coordinate (m)');
    legend('Perimeter between Wall and Channel - Uwc','Perimeter between Wall and Solid - Uws','Perimeter between Solid and Channel - Usc');
    figure; plot(A_T,Z,'b','LineWidth',3); hold on; plot(A_D,Z,'r','LineWidth',3); grid on;
    xlabel('Angular components (degrees)'); ylabel('y - vertical coordinate (m)'); legend('ThetaC','DeltaC');
end
if isempty(find(draw == 4, 1)) == 0 % Filled surface cylindrical normal pressure plot - 3D
    figure4 = figure('PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
    axes4 = axes('CameraPosition',[-59.84 -54.83 -17.81],'CameraUpVector',[-0.2156 -0.1975 0.9563],'OuterPosition',[0 0 1 1],'Parent',figure4);
    disp('Plotting 3D surface.'); surf4 = surf(srf.X,srf.Y,srf.Z,srf.P*Ch/1e3,'Parent',axes4); 
    grid(axes4,'on'); xlabel(axes4,'x'); ylabel(axes4,'y'); zlabel(axes4,'z'); axis equal;
    c = colorbar('Location','SouthOutside'); xlabel(c,'Normal pressures on silo wall, qn [kPa]'); shading interp; alpha(0.5);
end
if isempty(find(draw == 5, 1)) == 0; Font = 18; % Straight axis surface normal pressure plot - 3D
    figure5 = figure('Color',[0.9725 0.9725 0.9725],'Colormap',[0 0 0],'PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
    axes5 = axes('FontSize',Font,'Parent',figure5,'Position',[0.13 0.125 0.775 0.8783],...
        'ZTickLabel',{'1.0','0.8','0.6','0.4','0.2','0.0'}); 
    [dif, pos] = min(abs(srf.T(1,:) - Tmax));
    surf(srf.T(:,1:pos),srf.P(:,1:pos)*Ch/1e3,srf.Z(:,1:pos)/H,'FaceColor',[1 1 1],...
        'LineStyle',':','LineWidth',0.3,'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9608 0.9216 0.9216],...
        'MarkerSize',0.5,'FaceLighting','none','EdgeLighting','flat','Parent',axes5);
    %surf(srf.T,srf.P/1e3,srf.Z,'FaceColor',[1 1 1],'LineStyle',':','LineWidth',0.3, 'MarkerEdgeColor',[0 0 0],'MarkerFaceColor',[0.9608 0.9216 0.9216],'MarkerSize',0.5,'FaceLighting','none','EdgeLighting','flat');
    xlabel('Circumferential spread from 0 position (degrees)'); shading flat; view(210,30);
    ylabel('Normal pressure (kPa)'); zlabel('Normalised height, y/H'); 
end
if isempty(find(draw == 6, 1)) == 0 % 3D plot of the actualy geometry of channels and silo
    figure6 = figure('Colormap',[0 0 0],'PaperPosition',[0.6345 6.345 20.3 15.23],'PaperSize',[20.98 29.68]);
    axes6 = axes('CameraPosition',[25.29 -35.46 36.85],'CameraUpVector',[-0.3413 0.4785 0.809],...
        'XGrid','on','XTick',[-1 -0.5 0 0.5 1],'YGrid','on','YTick',[-1 -0.5 0 0.5 1],'ZTick',[0 H/R],'Parent',figure6);
    axis(axes6,[-1 1 -1 1 0 H/R]); view(axes6,[-37.5 30]); hold(axes6,'all');
    % 3D surface plots of the channel and silo
    surf(ch.X/R,ch.Y/R,ch.Z/R,'Parent',axes6); shading flat; hold on; 
    surf(st.X/R,st.Y/R,st.Z/R,'Parent',axes6); shading flat; colormap([0 0 0]); alpha(0.25); axis equal; 
    % Top surface circle
    [x1,y1,z1]=cylinder(R,100); plot(x1(1,:)/R,y1(1,:)/R,'k'); plot3(x1(1,:)/R,y1(1,:)/R,(H/R)*ones(length(x1(1,:))),'k');
     % Outlet circle
    [x2,y2,z2]=cylinder(rcm,100); plot(ec/R+x2(1,:)/R,y2(1,:)/R,'k');   
    % Symmetry plane
    plot3([R R -R -R R]/R,[0 0 0 0 0],[0 H H 0 0]/R,'k--'); 
    % Axis through channel centre
    plot3([ec ec]/R,[0 0]/R,[0 H]/R,'k-.');    
    % Region boundary citcle
    if z23 > H; plot3(x1(1,:)/R,y1(1,:)/R,(H-z23)*ones(length(x1(1,:)))/R,'k--'); end
    if z12 > H; plot3(x1(1,:)/R,y1(1,:)/R,(H-z12)*ones(length(x1(1,:)))/R,'k--'); end
    % Effective transition lines
    plot3(lc.x/R,lc.y/R,lc.z/R,'k','LineWidth',1.5); plot3(lc.x/R,-lc.y/R,lc.z/R,'k','LineWidth',1.5);
end


if scribe == 1; fprintf(fid,'%s\n',' '); fprintf(fid,'%s\n','#END'); fclose(fid); end
if scribe == 1; W = what; W.path(length(W.path)+1) = '\'; 
    W.path(length(W.path)+1:length(W.path)+length(File)) = File; D = dir(W.path);
    if D.bytes/1024 < 1e3; bytes = D.bytes/1024; ss = 'kB';
    else bytes = D.bytes/1024e3; ss = 'MB';
    end
    disp(['File written: ',File,' (',num2str(bytes),' ',ss,')']); 
end
% Thank you.
