function [] = AJSABQ_draw_2dplanar(Vals, S)
% Function to draw the planar ecentric discharge flow channel (if  applicable)
% AJSABQ Input File Generator Version 2.0
% Last modified by Dr Adam Jan Sadowski - 15/06/22 - 12:47 (previously 09/01/10 - 21:08)

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

kc = Vals.kc; ec = Vals.ec; Psi = Vals.Psi; THc = Vals.THc; Ac = Vals.Ac; 

% Plot outer circle
Circ_out.x = []; Circ_out.y = []; I1 = 0;
for I = 0:100 
    t = I*(2*pi/100);
    I1 = I1 + 1; Circ_out.x(I1) = cos(t); Circ_out.y(I1) = sin(t);    
end
% Plot innter circle
Circ_in.x = []; Circ_in.y = []; I2 = 0;
for I = 0:100 
    t = I*(2*pi/100);
    if t*180/pi <= Psi || t*180/pi >= (360-Psi); continue
    else
        I2 = I2 + 1; Circ_in.x(I2) = kc*cos(t)+ec; Circ_in.y(I2) = kc*sin(t);
    end
end

figure; plot(Circ_out.x,Circ_out.y,'k','LineWidth',3); axis square; grid on; hold on;
plot(Circ_in.x,Circ_in.y,'r','LineWidth',3); xlabel('X coordinate'); ylabel('Y coordinate');
line([-1 1],[0 0],'Color',[0 0.498 0],'LineWidth',2);
line([ec ec],[-1 1],'Color',[0 0.498 0],'LineWidth',2);
line([0 0],[-1 1],'Color',[0 0.498 0],'LineWidth',2);
line([0 cos(THc*pi/180)],[0 sin(THc*pi/180)],'Color',[0 0 1],'LineStyle',':','LineWidth',2); 
line([ec cos(THc*pi/180)],[0 sin(THc*pi/180)],'Color',[0 0 1],'LineStyle',':','LineWidth',2);
title(['Step ',num2str(S),'; Eccentric discharge planar geometry']);
annotation('textbox','Position',[0.03052 0.8158 0.0152 0.06084],'LineStyle','none','FitHeightToText','on','String',...
    {['kc = rc/R = ',num2str(kc)]});
annotation('textbox','Position',[0.03052 0.7158 0.0152 0.06084],'LineStyle','none','FitHeightToText','on','String',...
    {['ec/R = ',num2str(ec)]});
annotation('textbox','Position',[0.03052 0.6158 0.0152 0.06084],'LineStyle','none','FitHeightToText','on','String',...
    {['Ac/A = ',num2str(Ac*100),' %']});
annotation('textbox','Position',[0.03052 0.5158 0.0152 0.06084],'LineStyle','none','FitHeightToText','on','String',...
    {['ThetaC = ',num2str(THc),' degrees']});
annotation('textbox','Position',[0.03052 0.4158 0.0052 0.06084],'LineStyle','none','FitHeightToText','on','String',...
    {['PsiC = ',num2str(Psi),' degrees']});