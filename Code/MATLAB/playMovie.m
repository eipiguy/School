function [ fig ] = playMovie( mov )
%% playMovie: Summary of this function goes here
%   Detailed explanation goes here
% input:
% output:
%##########################################################################
%% Pseudo-Code:
%   ####
%   Figure Setup:
%   ====
%   Frame Render:
%   ====
%   Display:
%   ####
%##########################################################################
%% Figure Setup:

%scrsz = get(groot,'ScreenSize');
%set(gcf,'Position',[scrsz(3)/6, scrsz(4)/6, (scrsz(3))/3, (scrsz(4))/3]);

%==========================================================================
%% Frame Render:

function redraw(frame)
    imshow(mov(frame).cdata);
end

%==========================================================================
%% Display:

vfig = videofig(length(mov),@redraw);
redraw(1);

%##########################################################################
end

