clear all;
close all;

% set up path to .mat files
path ='./WindFields/';


% searching for files ending on '.mat'
filePattern = fullfile(path, '*.mat');
AllFiles = dir(filePattern);

for ind = 1:size(AllFiles,1)
    file = [AllFiles(ind).folder,'/',AllFiles(ind).name];
    disp (['converting file : ',file]);
    [fpath,fname,fext ] = fileparts(file);

    % It is assumed that the .mat file has a structure of [realizations, NumComp, ny, nz, t]
    load(file);

    n_realizations=size(u_gauss,1);
    n_components=size(u_gauss,2);
    ny=size(u_gauss,3);
    nz=size(u_gauss,4);
    nx=size(u_gauss,5);


    disp(['Number of realizations= ' ,num2str(n_realizations),...
        '  Components of wind=',num2str(n_components),...
        ' Size u=[',num2str(ny),'x',num2str(nz),'x',num2str(nx),']']);

    % !!!!!!!!!!!!!! Change name depending of 'gauss', 'temporal', 'Spatiotemporal'
     u_All=u_gauss;

    % ------------- Saving parameters of the field ------
    % Creating struct for parameters
    u_All_param=struct();
    % Sving parameters in ()
    for i_param=1:length(u_param)

        % Extracting name field from 'u_param_names' and value from 'u_param'
        u_All_param.(strtrim(u_param_names(i_param,:)))=u_param{i_param};
    end

    % 
    for realization_i=1:size(u_All,1) % In case of more than one realization
        
        %------------- Creating 4D array with wind field-----
        % Creating matrix 'u' [(time, 3D-windcomp, y, z)]
        u=zeros(nx,n_components,ny,nz);

        % Transforming [n_components, ny, nz, nx] -> [nx,n_components, ny, nz] 
        u=permute(squeeze(u_All(realization_i,:,:,:,:)),[4 1 2 3]);


        %-------- Formating to .wnd file
        %FileName: Name of output .wnd-file (extension will be added)
        % - zOffset: Reference height (m) = Z(1) + GridHeight / 2.0
        % - z0 = .03; Roughness length (m)
        % - z1 = Lowest location of the grid
        % - SummVars: 6 variables from the summary file {Clockwise, zHub, UBAR, TI_u, TI_v, TI_w}
    % TI_u, TI_v... Should be in % (i.e 10 and NOT 0.1)
        
        FileName=strcat(fpath,'/',fname,'_',num2str(realization_i));
        
        dy=u_All_param.y(2)-u_All_param.y(1);
        dz=u_All_param.z(2)-u_All_param.z(1);
        dt=u_All_param.T/nx;
       
        z0=0;
        z1=0;
        zHub=u_All_param.z(u_All_param.N_hub) - u_All_param.z(1);
        %zOffset=zHub;
        zOffset=z1 + ((nz-1)*dy) / 2.0;
        
        % Creating 'summVars' : Be careful with Turbulence Intensity!!
        % for example: if TI_v=0.8TI_u, then 80*(u_All_param.sigma/u_All_param.V_hub)

        TI_u=u_All_param.sigma/u_All_param.V_hub;

        SummVars=[1, zOffset, u_All_param.V_hub,...
            100*TI_u, 100*TI_u, 100*TI_u];

        % SummVars=[1, zHub, u_All_param.V_hub,...
        %     100*TI_u, 100*TI_u, 100*TI_u];

        % ------- Writing .wnd file-------- 
        writeBLgrid(FileName, u, dy, dz, dt, zOffset, z0, SummVars, z1)

        % ------- Writing .sum file-------- 
        writeSUM(FileName,SummVars(3),SummVars(4),SummVars(5),SummVars(6), zOffset, z1)
    end
end
