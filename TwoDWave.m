clear; clc;

myplot(20,0.5,1);

function myplot(nr,t, incr)
for i=0:incr:nr
    %
    %read the file simulation number <i>
    fn = sprintf('.u.dat.%03d',i);
    fp = fopen(fn,'r');
    d = fscanf(fp, '%f %f %f', [2 inf]);
    fclose(fp);
    %
    %plot the result;
    surf(d(1,:), d(2,:), d(3,:));
    axis([0 1 -0.1 0.1]);
    tittel = sprintf('Time t = %.3f', (t * i)/nr);
    title(tittel);
    
    %force drawing explicitly and wait 0.2 seconds
    drawnow;pause(0.5);
end
end