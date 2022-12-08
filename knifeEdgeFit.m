%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Knife edge fitting program
%
% Version: 1.2
% Date: 06/21/19
%
% Program performs fitting of knife edge data to determine the beam
% waist and center position of Gaussian beam. The input is a text file
% with the first column being the position and second column being the
% measured signal. It outputs the data that was input into the program
% the fitted data in one file, and the beam waist and central position in
% another file.
%
%     Copyright (C) 2017  Manuel Ferdinandus (mrf@knights.ucf.edu)
% 
%     This program is free software: you can redistribute it and/or modify
%     it under the terms of the GNU General Public License as published by
%     the Free Software Foundation, either version 3 of the License, or
%     (at your option) any later version.
% 
%     This program is distributed in the hope that it will be useful,
%     but WITHOUT ANY WARRANTY; without even the implied warranty of
%     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%     GNU General Public License for more details.
% 
%     You should have received a copy of the GNU General Public License
%     along with this program.  If not, see <http://www.gnu.org/licenses/>.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% import the knife edge data
current_directory = pwd; % the current directory
save_filepath = [current_directory,'\lastpath.mat'];

% look for the lastpath.mat file
if exist(save_filepath,'file') == 2
    
    last_path = open(save_filepath);
    [filename,pathname] = uigetfile(strcat(last_path.pathname,'*.txt'));

else
    
    last_path = pwd;    
    [filename,pathname] = uigetfile(strcat(last_path,'\*.txt'));
    
end

% perform the fitting only if the filename is valid
if ~strcmp(num2str(filename),'0') && ~strcmp(num2str(pathname),'0')
    
    save(save_filepath,'pathname');
    knifeEdge = importdata(strcat(pathname,filename),'\t',1); % Import the file

    % select the columns for the Z position and the voltage
    Z_column = 1;
    V_column = 2;
    startRow = 1; % the row at which the good data starts

    % extract the Z-position and the voltage data
    Z = knifeEdge.data(startRow:end,Z_column);
    V = knifeEdge.data(startRow:end,V_column);

    % sort the rows from smallest to largest
    [Z,sortIndex] = sortrows(Z);
    V = V(sortIndex);

    % normalized the voltage signal
    minSub = (V - min(V));
    T = minSub/max(minSub);

    % determine whether increasing or decreasing
    if T(1) > T(end)

        direction = -1;

    else

        direction = +1;

    end

    % find the center of the Gaussian
    ZcenterInd = find(T > 0.5);

    if direction == -1

        Zcenter = Z(ZcenterInd(end));

        % estimate the beam waist
        ZleftInd = find(T > 0.9);
        Zleft = Z(ZleftInd(end));
        ZrightInd = find(T < 0.1);
        Zright = Z(ZrightInd(1));
        w0 = (Zright - Zleft)/2;

    else

        Zcenter = Z(ZcenterInd(1));

        % estimate the beam waist
        ZleftInd = find(T < 0.1);
        Zleft = Z(ZleftInd(end));
        ZrightInd = find(T > 0.9);
        Zright = Z(ZrightInd(1));
        w0 = abs((Zright - Zleft)/2);

    end

    % fit the curve
    modelFun = @(p,x) (p(1) + direction*p(2)*(erf((sqrt(2)/p(3))*(x-p(4))))); % fitting function
    startingVals = [0.5,0.5,w0,Zcenter]; % inititalize the parameters
    [coefEsts,R,J,CovB,MSE] = nlinfit(Z,T,modelFun,startingVals);
    fitline = modelFun(coefEsts,Z);

    % rescale the data so that the fit line is at 1
    T_scaled = T./max(fitline);
    fitline_scaled = fitline./max(fitline);

    % display the raw data
    fit_plot = figure('Name','Knife edge fit');
    hold on
    plot(Z,T_scaled,'ob');
    plot(Z,fitline_scaled,'-r');
    title(['Knife edge fit of ',filename]);
    xlabel('Position (mm)');
    ylabel('Normalized signal');
    axis tight;
    box on;
    hold off;

    % extract the fitted parameters
    w0fit = coefEsts(3);
    Zcenterfit = coefEsts(4);
    
    stErr = sqrt(diag(CovB));
    w0err = stErr(3);
    x0err = stErr(4);
    RMSE = sqrt(MSE);

    % display messages
    disp(['For file ',filename,'.']);
    disp(['Spot size w (HW1/e^2M) is: ',num2str(w0fit),' mm.']);
    disp(['Standard error in w is: ',num2str(w0err),' mm (',num2str(100*w0err/w0fit),'%).']);
    disp(['Center position x_0 is: ',num2str(Zcenterfit),' mm.']);
    disp(['Standard error in x_0 is: ',num2str(x0err),' mm (',num2str(100*x0err/Zcenterfit),'%).']);
    disp(['RMSE is: ',num2str(RMSE),'.']);

    % export the data if requested
    dlgText = 'Export knife edge fit';
    doExport = questdlg([dlgText,'?'],'Export fit');
    if strcmp(doExport,'Yes')

        % get the export filename
        splitInd = regexp(filename,'.txt');
        exportFileName = strcat(filename(1:splitInd-1),'_knife_edge.txt');
        exportFilePath = [pathname,exportFileName];
        [exportFileName,exportPath] = uiputfile(exportFilePath,dlgText);
        exportFilePath = [exportPath,exportFileName];

        % generate the export data
        export_data = cat(2,Z,T_scaled,fitline_scaled);

        % write the data to the file
        export_fid = fopen(exportFilePath,'w+');
        fprintf(export_fid, 'Position (mm),Trans (data),Trans (fit)');
        dlmwrite(exportFilePath,export_data,'newline','pc','delimiter',',','precision',8,'-append','roffset',1);

        fclose(export_fid); % Close log file.

        % write the log file
        logFileName = strcat(filename(1:splitInd-1),'_log.txt');
        logFilePath = [exportPath,logFileName];

        log_fid = fopen(logFilePath,'w+');

        fprintf(log_fid,'%s\r\n\r\n',[datestr(now),': Begin knife edge fit program: ',mfilename('fullpath'),'.']);
        fprintf(log_fid,'%s\r\n\r\n',['Import data file: ',[pathname,filename],'.']);

        fprintf(log_fid,'%s\r\n','Fitted beam parameters:');
        fprintf(log_fid,'%s\r\n',['    w (HW1/e^2M):        ',num2str(w0fit),' (mm)']);
        fprintf(log_fid,'%s\r\n',['    w_err:               ',num2str(w0err),' (mm), ',num2str(100*w0err/w0fit),'%']);
        fprintf(log_fid,'%s\r\n',['    x0:                  ',num2str(Zcenterfit),' (mm)']);
        fprintf(log_fid,'%s\r\n',['    x0_err:              ',num2str(x0err),' (mm), ',num2str(100*x0err/Zcenterfit),'%']);
        fprintf(log_fid,'%s\r\n\r\n',['    RMSE:                ',num2str(RMSE)]);
        
        % save a picture of the plot
        splitInd = regexp(exportFilePath,'.txt');
        pictureFileName = strcat(exportFilePath(1:splitInd - 1),'.jpg');
        saveas(fit_plot,pictureFileName,'jpg')
        
        pictureFileName = strcat(exportFilePath(1:splitInd - 1),'.emf');
        saveas(fit_plot,pictureFileName,'emf')

        fprintf(log_fid,'%s\r\n\r\n',['Fit exported to: ',exportFilePath,'.']);
        
        fclose(export_fid);

    end;

end