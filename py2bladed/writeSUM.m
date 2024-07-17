%%-------------- This function was developed by Stefan Kapp (NREL) in 2009 https://forums.nrel.gov/t/binary-wnd-specification/180/5
% It has been modified by Daniela Moreno to match new
% requirements/definitions of FASTV8, OpenFast....
% --------------------------------------

function writeSUM(fileSum,u0_HH,TI_u,TI_v,TI_w, zOffset, z1)
%
% writes .sum files
%

% inputs - FileSum: Name of output .sum-file
%        - u0_HH: Mean wind speed at hub height
%        - TI_u: Turbulence Intensity of u-component at hub height in %
%        - TI_v: Turbulence Intensity of v-component at hub height in %
%        - TI_w: Turbulence Intensity of w-component at hub height in %

%-----------------------------------------
% OPEN FILE
%-----------------------------------------
sumfile=fopen(strcat(fileSum,'.sum'),'w');

        if sumfile>=0
            fprintf(sumfile,'%s %s\n', 'True', 'CLOCKWISE');
            % 'HUB HEIGHT' IS NOT REALLY the height of the hub but the
            % center of the grid in vertical direction or zOffset. Z(1) + GridHeight / 2.0
            fprintf(sumfile,'%2.2f %s\n\n',zOffset, 'HUB HEIGHT');
            fprintf(sumfile,'%2.5f %s\n',u0_HH, 'UBAR');
            fprintf(sumfile,'%2.5f%% %s\n',TI_u, 'TI(u)');
            fprintf(sumfile,'%2.5f%% %s\n',TI_v, 'TI(v)');
            fprintf(sumfile,'%2.5f%% %s\n',TI_w, 'TI(w)');
            

            fprintf(sumfile,'%2.2f %s\n\n',z1, 'GRID BASE');
            fclose(sumfile);
        else
            error('Could not open file for writing');
        end

end 