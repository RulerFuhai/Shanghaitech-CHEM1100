% MATLAB script to generate an AVI of a 1D traveling wave
% After Ralph Muehleisen, 2001
% 
%
%
%

% suppress warning
warning('off', 'MATLAB:avifile:indeo5NotFound');

% Input file name
[file,path] = uiputfile('standing_wave.avi');

% Delete existing file
filename = [path,file];
if exist(filename)
    delete(filename)
end
    
% Parameters
load('param.mat');
lambda = param.lambda;
A = param.A;

% Period
T = 1/ param.nu;

% Length of moive in seconds
len = 5;

% Position axis
x = (0 : 0.05 : 5) * lambda;

% create figure and size axes so we can get the geometry
% size for the movie
figure

axis([0, max(x),-A, A]);

set(gca,'nextplot','replacechildren');
% set up movie parameters
% use highest quality and don't compress as it introduces artifacts
% when we convert to a gif file we can compress that with less problems
% set frame rate to 7 frames per second so we can actually see what is happening
mov = VideoWriter(filename);
mov.FrameRate =7;
mov.Quality = 80;

% Number of movie frames
Nframes = mov.FrameRate * len;

open(mov);
for k = 1 : Nframes
    
    % time
    t = k / mov.FrameRate;
    
    % Traveling wave
    y = A * cos(2 * pi * (x / lambda)) * cos(2 * pi * (t / T));
	plot(x, y);

    title('Standing Wave')
    xlabel('Position')
    ylabel('Displacement');

    % actually get the movie frame and save it to the movie
    Frame = getframe(gcf);
    writeVideo(mov,Frame);
end
close(mov);

clearvars
