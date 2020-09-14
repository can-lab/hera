%----------------------------------------------------------------
%HERA: Heart Rate Analysis v0.07 beta April 2020
%
%(c) 2008-2020 Erno Hermans & Thijs vd Laar
%----------------------------------------------------------------

%v0.02: fixed check for wavelet manager
%       fixed use of nanmean (statistics toolbox): replaced by selfmade
%v0.03: fixed samplerate error on line 364 (read SR from file, not
%       defaults) 20110303 EJH
%v0.04: Added marker channel to visualize scanner triggers,
%       eg for Brainamp converted files. (20130310)
%v0.05: Changed behavior of peak correction. Now puts in and removes peaks
%       at the center of the figure, and does so much faster.
%v0.06: Added transparency to the red rectangles that indicate
%       the rejection period.
%v0.07: Added root median square of successive differences (rMedSSD) as a
%       potentially more robust measure of HRV. Updated "collect for output
%       to work also with incomplete data.
%v0.08: Replaced root median square of successive differences by root 
%       truncated mean of successive differences (rtMSSD), because the low sampling
%       frequency caused the median to jump.

function varargout = hera(varargin)
% HERA M-file for hera.fig
%      HERA, by itself, creates a new HERA or raises the existing
%      singleton*.
%
%      H = HERA returns the handle to a new HERA or the handle to
%      the existing singleton*.
%
%      HERA('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in HERA.M with the given input arguments.
%
%      HERA('Property','Value',...) creates a new HERA or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before hera_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to hera_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help hera

% Last Modified by GUIDE v2.5 02-Jun-2020 19:27:34

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @hera_OpeningFcn, ...
                   'gui_OutputFcn',  @hera_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before hera is made visible.
function hera_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to hera (see VARARGIN)

% Choose default command line output for hera
handles.output = hObject;




%------------------------------------------------------------------------
%Enter variables that need to be shared by subs and put them in "handles"
%------------------------------------------------------------------------

%Load Defaults
hera_defaults;
handles.hera_def=hera_def;                                  %these are default settings
handles.currentfile.matfile.settings = handles.hera_def;    %current settings
handles.currentfile.open = false;                           %currently no file open
handles.currentfile.saved = true;
clear hera_def;

%Check if complex wavelet (cmrl) has been added in wavemanager
cwavs = wavemngr('read');
cwavs = cwavs';
if isempty(strfind(cwavs(:)','cmrl'))
    %Add complex wavelet toolbox
    wavemngr('add', 'CMorlet', 'cmrl', 5, '', 'cmorwavf', [-3, 3]);
end


%Empty the matfile structure
handles = emptyMatfile(handles);
handles = drawgui(handles);

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes hera wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = hera_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



%----------------------------------------------------------------
%Open a file
%----------------------------------------------------------------




% --- Executes on button press in hera_openfile.
function hera_openfile_Callback(hObject, eventdata, handles)
% hObject    handle to hera_openfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if ~handles.currentfile.saved & handles.currentfile.open
    %Give a warning if current file has not been saved
    
    answ=questdlg('Current datafile has not been saved. Continue?');
    if isempty(strmatch('Yes',answ))
       
        %User cancelled, go back to gui
        return
        
    end
end

%Open a dialog box for opening a file
[pulsfilename,pulsdir,nocancel]=uigetfile({['*.',handles.hera_def.pulsefile_extension,';*.',upper(handles.hera_def.pulsefile_extension)],'Pulse files'},'Please specify the ".puls" file to work on.');

%If user pushes cancel, go back to gui
if ~nocancel
    return
end

%Empty possible current matfile structure
handles = emptyMatfile(handles);

%Load the pulse file: only read the first line of the file and convert to
%a vector
fid = fopen(fullfile(pulsdir,pulsfilename));
pulsfile = fread(fid);

%Find end of lines (ascii 13, 10)
eol = find((pulsfile(1:end-1)==13) & pulsfile(2:end)==10);
if numel(eol)>0
    eol = eol(1);
    pulsfile = pulsfile(1:eol-1);
end

%Now convert to numerical
handles.currentfile.pulsfile = str2num(char(pulsfile'));
handles.currentfile.pulsdir = pulsdir;
handles.currentfile.pulsfilename =pulsfilename;
handles.currentfile.open = true;

%Copy the pulsfile to the matfile (will later be edited)
handles.currentfile.matfile.rawpulsedata = handles.currentfile.pulsfile;

%Check if there already is a '.mat' file associated with this file.
%if so, load it.
extloc = strfind(handles.currentfile.pulsfilename,['.',handles.hera_def.pulsefile_extension]);
handles.currentfile.matfilename = [handles.currentfile.pulsfilename(1:extloc(length(extloc))),'mat'];
if exist(fullfile(handles.currentfile.pulsdir,handles.currentfile.matfilename))
    
    %Load the matfile and put it in the handles structure
    matfile = load(fullfile(handles.currentfile.pulsdir,handles.currentfile.matfilename));
    handles.currentfile.matfile = matfile.matfile;
else
    %Perform the first simple calculations on the pulse data
    
    %in case a pre-peak setting is present, use this to make an
    %ibichannel based on existing peak detection
    
    if handles.hera_def.prepeak > 0
        
        prepeakbool = handles.currentfile.matfile.rawpulsedata==handles.hera_def.prepeak;
        prepeaklocs = find(prepeakbool);
        
        if ~isempty(prepeaklocs)
            %Correct for Siemens-inserted sample points per peak
            prepeaklocs = prepeaklocs - [0:length(prepeaklocs)-1];

            %Correct for peak detection delay in Siemens scanner data,
            %according to setting in defaults
            prepeaklocs = prepeaklocs - handles.hera_def.prepeakdelay;

            %Store the pre-peak locations for later use
            handles.currentfile.matfile.prepeaklocs = prepeaklocs;
        
            %Now delete all peaks from the pulse data
            handles.currentfile.matfile.rawpulsedata = ...
                handles.currentfile.matfile.rawpulsedata(~prepeakbool);
        
        
        %Loop over all peaks, if any, and delete all the peaks from the raw
        %data
       % if ~isempty(prepeaklocs)
            
       %     for cpeak = 1:length(prepeaklocs)
                
                %Interpolate the datapoint
       %         if prepeaklocs(cpeak) >1 & prepeaklocs(cpeak)<length(handles.currentfile.matfile.rawpulsedata)
                    
                    %Simply interpolate from neighbors
       %             handles.currentfile.matfile.rawpulsedata(prepeaklocs(cpeak)) = mean ([...
       %                 handles.currentfile.matfile.rawpulsedata(prepeaklocs(cpeak)-1), ...
       %                 handles.currentfile.matfile.rawpulsedata(prepeaklocs(cpeak)+1)]);
                    
       %         elseif prepeaklocs(cpeak) ==1
                    
                    %If it's the first datapoint, copy the second
       %             handles.currentfile.matfile.rawpulsedata(1) = handles.currentfile.matfile.rawpulsedata(2);

       %         elseif prepeaklocs(cpeak) == length(handles.currentfile.matfile.rawpulsedata)
                    %If it's the last datapoint, copy the previous
       %             handles.currentfile.matfile.rawpulsedata(length(handles.currentfile.matfile.rawpulsedata)) = ...
       %                 handles.currentfile.matfile.rawpulsedata(length(handles.currentfile.matfile.rawpulsedata)-1);
       %         end
                
       %     end
            
            
            %Now filter the data according to the filter settings, in
            %this case the default settings because no file is open

            handles.currentfile.matfile.bpf_pulsedata = ...
                hera_bpf(handles.currentfile.matfile.rawpulsedata, ...
                    handles.currentfile.matfile.settings.samplerate, ...
                    handles.currentfile.matfile.settings.hpd, ...
                    handles.currentfile.matfile.settings.lpd);

            %Define time-axis
            handles.currentfile.matfile.sampleTime = ...
                (1/handles.currentfile.matfile.settings.samplerate) * ...
                (0:length(handles.currentfile.matfile.bpf_pulsedata) - 1);
            
            %Determine the time (in sec) at which beats were detected
            handles.currentfile.matfile.prepeakTimes = ...
                handles.currentfile.matfile.sampleTime(prepeaklocs);
            
            %Call other function to calculate IBI timeseries and channels
            handles = CalcIBIs(handles);
            
            %Calculate rMSSD and BPM measures, use BPM as pre-setting for
            %wavelet peak tracing.
            handles = CalcPreOutMeasures(handles);
        end
    end
end


handles = drawgui(handles);

guidata(hObject, handles);




%----------------------------------------------------------------
%Save a file
%----------------------------------------------------------------




% --- Executes on button press in hera_save.
function hera_save_Callback(hObject, eventdata, handles)
% hObject    handle to hera_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Check if there is a datafile open
if handles.currentfile.open
    
    %Save everything in the matfile structure
    matfile = handles.currentfile.matfile;
    save(fullfile(handles.currentfile.pulsdir,handles.currentfile.matfilename),'matfile');
    handles.currentfile.saved = true;
else
    msgbox('No file has been opened. Nothing to save.');
    return
end

guidata(hObject, handles);    
    
    
%----------------------------------------------------------------
%Draw the GUI
%----------------------------------------------------------------
function handles = drawgui(handles)

%Calculate BPM/rMSSD measures
handles = CalcIBIs(handles);
handles = CalcPreOutMeasures(handles);

%Fill in all settings/results fields; these fields contain defaults
%when no file is open (or are empty)
set(handles.hera_lowpassfilter, 'String',handles.currentfile.matfile.settings.lpd);
set(handles.hera_highpassfilter, 'String',handles.currentfile.matfile.settings.hpd);
set(handles.hera_avgHeartbeat, 'String',handles.currentfile.matfile.settings.avgHeartbeat);
set(handles.hera_ballmobilityrange, 'String',handles.currentfile.matfile.settings.ballmobilityrange);
set(handles.hera_startCutout, 'String',handles.currentfile.matfile.settings.startCutout);
set(handles.hera_endCutout, 'String',handles.currentfile.matfile.settings.endCutout);
set(handles.hera_lowerboundLF, 'String',handles.currentfile.matfile.settings.lowerboundLF);
set(handles.hera_upperboundLF, 'String',handles.currentfile.matfile.settings.upperboundLF);
set(handles.hera_lowerboundHF, 'String',handles.currentfile.matfile.settings.lowerboundHF);
set(handles.hera_upperboundHF, 'String',handles.currentfile.matfile.settings.upperboundHF);
set(handles.hera_preBPM, 'String',handles.currentfile.matfile.outmeasures.preBPM);
set(handles.hera_prerMSSD, 'String',handles.currentfile.matfile.outmeasures.rMSSD);
set(handles.hera_prertMSSD, 'String',handles.currentfile.matfile.outmeasures.rtMSSD);
set(handles.hera_postBPM, 'String',handles.currentfile.matfile.outmeasures.postBPM);
set(handles.hera_LFHFrat, 'String',handles.currentfile.matfile.outmeasures.LFHFratio);
set(handles.hera_LFpower, 'String',handles.currentfile.matfile.outmeasures.LFPower);
set(handles.hera_HFpower, 'String',handles.currentfile.matfile.outmeasures.HFPower);



%Ylabel X popsition
Ylx = -15; 

%Reset progress bar
cposition = get(handles.hera_redprogressbar,'Position');
cposition(3) = 0.001;
set(handles.hera_redprogressbar,'Position',cposition);

%Reset the progress bar text
set(handles.hera_progressbartext,'String','');


%Draw the raw data
if handles.currentfile.open

    %Set the filename
    set(handles.figure1,'Name',['HeRA - ',handles.currentfile.pulsfilename]);
    
    %Select upper graph
    set(gcf,'CurrentAxes',handles.axes1);
    
    %If there's a pre-IBIchannel, plot in upper graph
    if ~isempty(handles.currentfile.matfile.preibichannel)
        
        %Rescale the pulsfile to fit in the same graph as the IBI channel
        minibi = min(handles.currentfile.matfile.preibichannel);
        maxibi = max(handles.currentfile.matfile.preibichannel);
        
        %20130316: set max to 2000 if artefacts are in the data
        if maxibi>2000
            maxibi=2000;
        end

        %Rescale the band pass filtered pulsefile
        rspulsfile = handles.currentfile.matfile.bpf_pulsedata;
        rspulsfile = rspulsfile.*((maxibi-minibi)/(max(rspulsfile) - min(rspulsfile)));
        rspulsfile = rspulsfile - min(rspulsfile);
        rspulsfile = rspulsfile + minibi;
        xaxis = ((1:length(rspulsfile))-1)./handles.currentfile.matfile.settings.samplerate; %20110303 changed to samplerate read from the file
        
        %Keep the x axis (20130316)
        handles.xaxis = xaxis;
        
        cYLim = [minibi-((maxibi-minibi)*.05),maxibi+((maxibi-minibi)*.05)];
        
        %Store the X limits so zooming can be checked later
        handles.currentfile.matfile.prepulsezoom = get(handles.axes1,'XLim');
        
        %20130316: should this be done later?
        plot(xaxis,rspulsfile,'k','linewidth',1)
        
        
        %Align Y axis labels
      %  ylpos = get(get(handles.axes1,'YLabel'),'Position');
     %   ylpos(1) = Ylx;
     %   set(get(handles.axes1,'YLabel'),'Position',ylpos);

        hold on
        %Draw in detected peaks
        handles.plothandle_peaks=[];
        for peak = 1:length(handles.currentfile.matfile.prepeakTimes)
            %Draw in vertical lines at peaks
            xPos = handles.currentfile.matfile.prepeakTimes(peak);
            %20130316: Keep handles to peaks in the figure
            handles.plothandle_peaks(peak) = ...
                plot([xPos,xPos],cYLim,'Color',[.8,.8,.8],'linewidth',.5);
        end

        %Draw in marker channel if present (added 20130310)
        if isfield(handles.currentfile.matfile,'markerlocs')
            markerYLim = [cYLim(1),cYLim(1)+((cYLim(2)-cYLim(1))/40)];
            for cmarkerloc=1:numel(handles.currentfile.matfile.markerlocs)
                xPos = handles.currentfile.matfile.markerlocs(cmarkerloc)/...
                    handles.currentfile.matfile.settings.samplerate;
                plot([xPos,xPos],markerYLim,'Color',[.5,.5,1],'linewidth',4)
            end
        end
        
        handles.plothandle_rspulsfile = plot(xaxis,rspulsfile,'k','linewidth',1);
        
        %Added 20131316: keep a handle to the preibichannel trace in the
        %plot so that it can be redrawn more easily
        handles.plothandle_preibichannel = plot(xaxis,handles.currentfile.matfile.preibichannel,'r','linewidth',1);
        
        set(get(gcf,'CurrentAxes'),'YLim',cYLim);
        
%         %added 20130316: return to previous zoom if available (EJH)
%         if isfield(handles,'keepzoom') && numel(handles.keepzoom)>0
%             set(get(gcf,'CurrentAxes'),'XLim',handles.keepzoom(1,:));
%             set(get(gcf,'CurrentAxes'),'YLim',handles.keepzoom(2,:));
%             handles.keepzoom=[];
%         else
%             set(get(gcf,'CurrentAxes'),'YLim',cYLim);
%         end    
        %Draw in rejected areas
        handles.plothandle_rejected=[];
        for rej = 1:length(handles.currentfile.matfile.prereject)
          
            %Draw in a rectangle to indicate rejection
            crej = handles.currentfile.matfile.prereject{rej};
            
            %20180201: added transparency to the red rectangles that indicate
            %the rejection period.
            handles.plothandle_rejected(rej) = ...
                patch(...
                    [crej(1) ,crej(2) ,crej(2) , crej(1)],...
                    [cYLim(1),cYLim(1),cYLim(2),cYLim(2)],...
                    [1,.7,.7]);
            set(handles.plothandle_rejected(rej),'FaceAlpha',.5);
            set(handles.plothandle_rejected(rej),'LineStyle','none');
    
            
            %old: using a rectangle which cannot be transparent:
%             handles.plothandle_rejected(rej) = ...
%                 rectangle('Position',...
%                 [crej(1),cYLim(1),crej(2)-crej(1),cYLim(2)-cYLim(1)]);
%             set(handles.plothandle_rejected(rej),'FaceColor',[1,.7,.7]);
%             set(handles.plothandle_rejected(rej),'LineStyle','none');
            
            %old: Draw in red vertical lines in rejected period
            %crej = handles.currentfile.matfile.prereject{rej};
            %for xPos = crej(1):.04:crej(2)
            %    plot([xPos,xPos],cYLim,'Color',[1,.7,.7])
            %end
        end
        hold off
        xlabel('Time(s)'), ylabel('InterBeat Interval(ms)')
    else
        
        %Rescale the raw pulsefile to 0-100
        rspulsfile = handles.currentfile.matfile.bpf_pulsedata;
        rspulsfile = rspulsfile.*(100/(max(rspulsfile) - min(rspulsfile)));
        rspulsfile = rspulsfile - min(rspulsfile);
        plot(handles.currentfile.matfile.sampleTime,rspulsfile,'k','linewidth',1.5)
    
    end
    
    %If there's a wavelet transform available, draw it
    if ~isempty(handles.currentfile.matfile.CWT)

        set(gcf,'CurrentAxes',handles.axes2);
        
        imagesc(handles.currentfile.matfile.sampleTime, ...             %Plot CWT against frequency and time
            handles.currentfile.matfile.scalFreq, ...
            handles.currentfile.matfile.CWT)
        xlabel('Time(s)'), ylabel('Frequency(Hz)'); %Proper labeling
        
        %Align Y axis labels
 %       ylpos = get(get(handles.axes2,'YLabel'),'Position');
 %       ylpos(1) = Ylx;
 %       set(get(handles.axes2,'YLabel'),'Position',ylpos);

    else
        set(gcf,'CurrentAxes',handles.axes2);
        cla
    end

    %If there's a top-trace available, draw it in
    if ~isempty(handles.currentfile.matfile.topTrace)
        hold on
        plot(handles.currentfile.matfile.sampleTime, ...
            handles.currentfile.matfile.topTrace, 'k');
        hold off
        xlabel('Time(s)'), ylabel('Frequency(Hz)'); %Proper labeling
    end
    
    
    
    if ~isempty(handles.currentfile.matfile.spectrumPower)
    

        startFreqIndex = round(handles.currentfile.matfile.lowerboundLFIndex - (handles.currentfile.matfile.upperboundHFIndex - handles.currentfile.matfile.lowerboundLFIndex)*.1);
        if startFreqIndex<1; startFreqIndex = 1; end
        endFreqIndex = round(handles.currentfile.matfile.upperboundHFIndex + (handles.currentfile.matfile.upperboundHFIndex - handles.currentfile.matfile.lowerboundLFIndex)*.1);
        if endFreqIndex>length(handles.currentfile.matfile.spectrumFrequencies)
            endFreqIndex = length(handles.currentfile.matfile.spectrumFrequencies);
        end

        sF = handles.currentfile.matfile.spectrumFrequencies(startFreqIndex:endFreqIndex);
        sP = handles.currentfile.matfile.spectrumPower(startFreqIndex:endFreqIndex);
        
        set(gcf,'CurrentAxes',handles.axes3);
        plot(sF, sP,'k','linewidth',2)
        
        
        %Align Y axis labels
  %      ylpos = get(get(handles.axes3,'YLabel'),'Position');
  %      ylpos(1) = -.049;
  %      set(get(handles.axes3,'YLabel'),'Position',ylpos);
        
        
        %Put in lines to see LF and HF boundaries
        lLF = handles.currentfile.matfile.settings.lowerboundLF;
        uLF = handles.currentfile.matfile.settings.upperboundLF;
        lHF = handles.currentfile.matfile.settings.lowerboundHF;
        uHF = handles.currentfile.matfile.settings.upperboundHF;
        
        cYLim = get(handles.axes3,'YLim');
        hold on
        for xPos = lLF:.001:uLF
            plot([xPos,xPos],cYLim,'Color',[.7,.7,1])
        end
        for xPos = lHF:.001:uHF
            plot([xPos,xPos],cYLim,'Color',[.7,1,.7])
        end
        plot(sF, sP,'k','linewidth',2)
        hold off
        xlabel('Frequency(Hz)'), ylabel('Power spectral density (W/Hz)');
                
    end
    
    
    
else %if no file open, clear all graphs
    
    %Select graphs and clear
    set(gcf,'CurrentAxes',handles.axes1);
    cla
    set(gcf,'CurrentAxes',handles.axes2);
    cla
    set(gcf,'CurrentAxes',handles.axes3);
    cla
    
end
    
%----------------------------------------------------------------
%Close a file
%----------------------------------------------------------------


% --- Executes on button press in hera_closefile.
function hera_closefile_Callback(hObject, eventdata, handles)
% hObject    handle to hera_closefile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check if there is a datafile open
if handles.currentfile.open
    
    %Check if file has been saved
    if ~handles.currentfile.saved
        %Give a warning if current file has not been saved
        answ=questdlg('Current datafile has not been saved. Continue?');
        if isempty(strmatch('Yes',answ))

            %User cancelled, go back to gui
            return

        end
    end
    
    %Close the file and redraw
    handles.currentfile.open = false;
    handles = emptyMatfile(handles);
    set(handles.figure1,'Name','HeRA');
    handles = drawgui(handles);
    guidata(hObject, handles);
end

%----------------------------------------------------------------
%Close a file and clear/delete the matfile
%----------------------------------------------------------------


% --- Executes on button press in hera_clearmatfile.
function hera_clearmatfile_Callback(hObject, eventdata, handles)
% hObject    handle to hera_clearmatfile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Only do something if a file is open
if handles.currentfile.open
    
    answ = questdlg('This will delete all changes, calculations, and the ".mat"-file of the currently opened file, and then close it. Are you sure?');
    if ~strmatch(answ,'Yes')
        return
    end
    
    %Close the file and delete the mat-file
    handles.currentfile.open = false;
    handles = emptyMatfile(handles);
    
    %Delete the matfile
    warning off;
    delete(handles.currentfile.matfilename);
    warning on;
    
    set(handles.figure1,'Name','HeRA');
    handles = drawgui(handles);
    guidata(hObject, handles);

end




%---------------------------------------------------
%Emptymatfile: create the empty matfile structure 

function handles = emptyMatfile(handles)

handles.currentfile.matfile.settings = handles.hera_def;    %go back to default settings
handles.currentfile.matfile.rawpulsedata = [];
handles.currentfile.matfile.preibitimeseries = [];
handles.currentfile.matfile.preibichannel = [];
handles.currentfile.matfile.prepeaklocs = [];
handles.currentfile.matfile.prepeakTimes = [];
handles.currentfile.matfile.prereject = [];
handles.currentfile.matfile.prepulsezoom = [];
handles.currentfile.matfile.postibichannel = [];
handles.currentfile.matfile.sampleTime = [];
handles.currentfile.matfile.scalFreq = [];
handles.currentfile.matfile.aBin = [];
handles.currentfile.matfile.CWT = [];
handles.currentfile.matfile.topTrace = [];
handles.currentfile.matfile.bpf_pulsedata = [];
handles.currentfile.matfile.spectrumFrequencies = [];
handles.currentfile.matfile.spectrumPower = [];
handles.currentfile.matfile.lowerboundLFIndex = [];
handles.currentfile.matfile.upperboundLFIndex = [];
handles.currentfile.matfile.lowerboundHFIndex = [];
handles.currentfile.matfile.upperboundHFIndex = [];
handles.currentfile.matfile.outmeasures.preBPM = [];
handles.currentfile.matfile.outmeasures.rMSSD = [];
handles.currentfile.matfile.outmeasures.rtMSSD = [];
handles.currentfile.matfile.outmeasures.LFPower = [];
handles.currentfile.matfile.outmeasures.HFPower = [];
handles.currentfile.matfile.outmeasures.LFHFratio = [];
handles.currentfile.matfile.outmeasures.postBPM = [];



%---------------------------------------------------



%------------------------------------------------------------------
%Change initial low pass filter
%------------------------------------------------------------------

function hera_lowpassfilter_Callback(hObject, eventdata, handles)
% hObject    handle to hera_lowpassfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_lowpassfilter as text
%        str2double(get(hObject,'String')) returns contents of hera_lowpassfilter as a double


if ~handles.currentfile.open   
    msgbox('No file has been opened. Cannot change this setting now.')
    set(handles.hera_lowpassfilter, 'String',handles.currentfile.matfile.settings.lpd);
    return
end
    

cset = get(handles.hera_lowpassfilter, 'String');

if ~isempty(str2num(cset))
    
    %Set the filter setting in the matfile
    handles.currentfile.matfile.settings.lpd = str2num(cset);
    
    %Refilter the raw pulse data    
    
    handles.currentfile.matfile.bpf_pulsedata = ...
    hera_bpf(handles.currentfile.matfile.rawpulsedata, ...
        handles.currentfile.matfile.settings.samplerate, ...
        handles.currentfile.matfile.settings.hpd, ...
        handles.currentfile.matfile.settings.lpd);

    
    %Redraw the GUI
    handles = drawgui(handles);
    
    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;
    
else
    msgbox('Invalid setting. Value set back to previous setting.')
    set(handles.hera_lowpassfilter, 'String',handles.currentfile.matfile.settings.lpd);
end

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function hera_lowpassfilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_lowpassfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%------------------------------------------------------------------
%Change initial high pass filter
%------------------------------------------------------------------



function hera_highpassfilter_Callback(hObject, eventdata, handles)
% hObject    handle to hera_highpassfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_highpassfilter as text
%        str2double(get(hObject,'String')) returns contents of hera_highpassfilter as a double

if ~handles.currentfile.open   
    msgbox('No file has been opened. Cannot change this setting now.')
    set(handles.hera_highpassfilter, 'String',handles.currentfile.matfile.settings.hpd);
    return
end

cset = get(handles.hera_highpassfilter, 'String');

if ~isempty(str2num(cset))
    
    %Set the filter setting in the matfile
    handles.currentfile.matfile.settings.hpd = str2num(cset);
    
    %Refilter the raw pulse data    
    
    handles.currentfile.matfile.bpf_pulsedata = ...
    hera_bpf(handles.currentfile.matfile.rawpulsedata, ...
        handles.currentfile.matfile.settings.samplerate, ...
        handles.currentfile.matfile.settings.hpd, ...
        handles.currentfile.matfile.settings.lpd);

    
    %Redraw the GUI
    handles = drawgui(handles);
    
    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;

else
    msgbox('Invalid setting. Value set back to previous setting.')
    set(handles.hera_highpassfilter, 'String',handles.currentfile.matfile.settings.hpd);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function hera_highpassfilter_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_highpassfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%------------------------------------------------------------------
%Change average heartbeat setting
%------------------------------------------------------------------

function hera_avgHeartbeat_Callback(hObject, eventdata, handles)
% hObject    handle to hera_avgHeartbeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_avgHeartbeat as text
%        str2double(get(hObject,'String')) returns contents of hera_avgHeartbeat as a double

if ~handles.currentfile.open   
    msgbox('No file has been opened. Cannot change this setting now.')
    set(handles.hera_avgHeartbeat, 'String',handles.currentfile.matfile.settings.avgHeartbeat);
    return
end

cset = get(handles.hera_avgHeartbeat, 'String');

if ~isempty(str2num(cset))
    
    %Set the filter setting in the matfile
    handles.currentfile.matfile.settings.avgHeartbeat = str2num(cset);
    
    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;

else
    msgbox('Invalid setting. Value set back to previous setting.')
    set(handles.hera_avgHeartbeat, 'String',handles.currentfile.matfile.settings.avgHeartbeat);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function hera_avgHeartbeat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_avgHeartbeat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%------------------------------------------------------------------
%Change ball mobility range setting
%------------------------------------------------------------------

function hera_ballmobilityrange_Callback(hObject, eventdata, handles)
% hObject    handle to hera_ballmobilityrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_ballmobilityrange as text
%        str2double(get(hObject,'String')) returns contents of hera_ballmobilityrange as a double

if ~handles.currentfile.open   
    msgbox('No file has been opened. Cannot change this setting now.')
    set(handles.hera_ballmobilityrange, 'String',handles.currentfile.matfile.settings.ballmobilityrange);
    return
end

cset = get(handles.hera_ballmobilityrange, 'String');

if ~isempty(str2num(cset))
    
    %Set the filter setting in the matfile
    handles.currentfile.matfile.settings.ballmobilityrange = str2num(cset);
    
    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;

else
    msgbox('Invalid setting. Value set back to previous setting.')
    set(handles.hera_ballmobilityrange, 'String',handles.currentfile.matfile.settings.ballmobilityrange);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function hera_ballmobilityrange_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_ballmobilityrange (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%------------------------------------------------------------------
%Change start cutout setting
%------------------------------------------------------------------

function hera_startCutout_Callback(hObject, eventdata, handles)
% hObject    handle to hera_startCutout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_startCutout as text
%        str2double(get(hObject,'String')) returns contents of hera_startCutout as a double

if ~handles.currentfile.open   
    msgbox('No file has been opened. Cannot change this setting now.')
    set(handles.hera_startCutout, 'String',handles.currentfile.matfile.settings.startCutout);
    return
end

cset = get(handles.hera_startCutout, 'String');

if ~isempty(str2num(cset))
    
    %Set the filter setting in the matfile
    handles.currentfile.matfile.settings.startCutout = str2num(cset);
    
    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;

else
    msgbox('Invalid setting. Value set back to previous setting.')
    set(handles.hera_startCutout, 'String',handles.currentfile.matfile.settings.startCutout);
end

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function hera_startCutout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_startCutout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


%------------------------------------------------------------------
%Change end cutout setting
%------------------------------------------------------------------

function hera_endCutout_Callback(hObject, eventdata, handles)
% hObject    handle to hera_endCutout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_endCutout as text
%        str2double(get(hObject,'String')) returns contents of hera_endCutout as a double

if ~handles.currentfile.open   
    msgbox('No file has been opened. Cannot change this setting now.')
    set(handles.hera_endCutout, 'String',handles.currentfile.matfile.settings.endCutout);
    return
end

cset = get(handles.hera_endCutout, 'String');

if ~isempty(str2num(cset))
    
    %Set the filter setting in the matfile
    handles.currentfile.matfile.settings.endCutout = str2num(cset);
    
    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;
    
else
    msgbox('Invalid setting. Value set back to previous setting.')
    set(handles.hera_endCutout, 'String',handles.currentfile.matfile.settings.endCutout);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function hera_endCutout_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_endCutout (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end










%------------------------------------------------------------------
%Change lowerbound LF setting
%------------------------------------------------------------------


function hera_lowerboundLF_Callback(hObject, eventdata, handles)
% hObject    handle to hera_lowerboundLF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_lowerboundLF as text
%        str2double(get(hObject,'String')) returns contents of hera_lowerboundLF as a double


if ~handles.currentfile.open   
    msgbox('No file has been opened. Cannot change this setting now.')
    set(handles.hera_lowerboundLF, 'String',handles.currentfile.matfile.settings.lowerboundLF);
    return
end

cset = get(handles.hera_lowerboundLF, 'String');

if ~isempty(str2num(cset)) & str2num(cset) > 0 & ...
        str2num(cset)< handles.currentfile.matfile.settings.upperboundLF
    
    %Set the filter setting in the matfile
    handles.currentfile.matfile.settings.lowerboundLF = str2num(cset);
    handles = CalcOutmeasures(handles);
    handles = drawgui(handles);
    
    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;
    
else
    msgbox('Invalid setting. Value set back to previous setting.')
    set(handles.hera_lowerboundLF, 'String',handles.currentfile.matfile.settings.lowerboundLF);
end

guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function hera_lowerboundLF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_lowerboundLF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------------------------------------------
%Change upperbound LF setting
%------------------------------------------------------------------


function hera_upperboundLF_Callback(hObject, eventdata, handles)
% hObject    handle to hera_upperboundLF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_upperboundLF as text
%        str2double(get(hObject,'String')) returns contents of hera_upperboundLF as a double


if ~handles.currentfile.open   
    msgbox('No file has been opened. Cannot change this setting now.')
    set(handles.hera_upperboundLF, 'String',handles.currentfile.matfile.settings.upperboundLF);
    return
end

cset = get(handles.hera_upperboundLF, 'String');

if ~isempty(str2num(cset)) & str2num(cset) > 0 & ...
        str2num(cset)> handles.currentfile.matfile.settings.lowerboundLF
    
    %Set the filter setting in the matfile
    handles.currentfile.matfile.settings.upperboundLF = str2num(cset);
    handles = CalcOutmeasures(handles);
    handles = drawgui(handles);
    
    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;
    
else
    msgbox('Invalid setting. Value set back to previous setting.')
    set(handles.hera_upperboundLF, 'String',handles.currentfile.matfile.settings.upperboundLF);
end

guidata(hObject, handles);





% --- Executes during object creation, after setting all properties.
function hera_upperboundLF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_upperboundLF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------------------------------------------
%Change lowerbound HF setting
%------------------------------------------------------------------


function hera_lowerboundHF_Callback(hObject, eventdata, handles)
% hObject    handle to hera_lowerboundHF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_lowerboundHF as text
%        str2double(get(hObject,'String')) returns contents of hera_lowerboundHF as a double


if ~handles.currentfile.open   
    msgbox('No file has been opened. Cannot change this setting now.')
    set(handles.hera_lowerboundHF, 'String',handles.currentfile.matfile.settings.lowerboundHF);
    return
end

cset = get(handles.hera_lowerboundHF, 'String');

if ~isempty(str2num(cset)) & str2num(cset) > 0 & ...
        str2num(cset)< handles.currentfile.matfile.settings.upperboundHF
    
    %Set the filter setting in the matfile
    handles.currentfile.matfile.settings.lowerboundHF = str2num(cset);
    handles = CalcOutmeasures(handles);
    handles = drawgui(handles);
    
    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;
    
else
    msgbox('Invalid setting. Value set back to previous setting.')
    set(handles.hera_lowerboundHF, 'String',handles.currentfile.matfile.settings.lowerboundHF);
end

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function hera_lowerboundHF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_lowerboundHF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end

%------------------------------------------------------------------
%Change upperbound HF setting
%------------------------------------------------------------------


function hera_upperboundHF_Callback(hObject, eventdata, handles)
% hObject    handle to hera_upperboundHF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_upperboundHF as text
%        str2double(get(hObject,'String')) returns contents of hera_upperboundHF as a double


if ~handles.currentfile.open   
    msgbox('No file has been opened. Cannot change this setting now.')
    set(handles.hera_upperboundHF, 'String',handles.currentfile.matfile.settings.upperboundHF);
    return
end

cset = get(handles.hera_upperboundHF, 'String');

if ~isempty(str2num(cset)) & str2num(cset) > 0 & ...
        str2num(cset)> handles.currentfile.matfile.settings.lowerboundHF
    
    %Set the filter setting in the matfile
    handles.currentfile.matfile.settings.upperboundHF = str2num(cset);
    handles = CalcOutmeasures(handles);
    handles = drawgui(handles);
    
    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;
    
else
    msgbox('Invalid setting. Value set back to previous setting.')
    set(handles.hera_upperboundHF, 'String',handles.currentfile.matfile.settings.upperboundHF);
end

guidata(hObject, handles);



% --- Executes during object creation, after setting all properties.
function hera_upperboundHF_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_upperboundHF (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end




%----------------------------------------------------------------
%Calculate CWT
%----------------------------------------------------------------


% --- Executes on button press in hera_calccwt.
function hera_calccwt_Callback(hObject, eventdata, handles)
% hObject    handle to hera_calccwt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.currentfile.open
    
    %Initialise variables and scale resolutions
    aRes = 0.00005;

    %Calculate corresponding freq resulutions to scales
    handles.currentfile.matfile.aBin = ...
        scal2frq(1/aRes, 'cmrl', 1/handles.currentfile.matfile.settings.samplerate);
    
    
    a = 1./(0.01:aRes:0.050); %Define scaling factors, may be adapted for resolution or range. Check corresponding frequency of scaling factor with scal2frq(). 0.5-2.5 Hz

    handles.currentfile.matfile.scalFreq = ...
        scal2frq(a, 'cmrl', 1/handles.currentfile.matfile.settings.samplerate); %Convert scaling factor to frequencies

    %Set the progress bar text to show what it's doing
    set(handles.hera_progressbartext,'String','Calculating CWT...');

    %Calculate the CWT and pass handles for progress bar access
    handles.currentfile.matfile.CWT = hera_cwt(handles.currentfile.matfile.bpf_pulsedata, a, 'cmrl',handles); %Complex Morlet waveform
    handles.currentfile.matfile.CWT = abs(handles.currentfile.matfile.CWT); %Length of imaginary vectors
    
    handles = drawgui(handles); 

    %Set the saved setting to "no" to avoid losing changes
    handles.currentfile.saved = false;

    guidata(hObject, handles);
end



%------------------------------------------------------------------
%Trace heartbeat on CWT
%------------------------------------------------------------------



% --- Executes on button press in hera_trace.
function hera_trace_Callback(hObject, eventdata, handles)
% hObject    handle to hera_trace (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Calculate the sample number to begin from and end at using the
%startCutout and endCutout settings
startFnd = find(handles.currentfile.matfile.sampleTime > ...
    handles.currentfile.matfile.settings.startCutout);
startM = startFnd(1);
endFnd = find(handles.currentfile.matfile.sampleTime > ...
    (handles.currentfile.matfile.sampleTime(end)- ...
    handles.currentfile.matfile.settings.endCutout));
endM = endFnd(1);



%Progress bar settings
pperc = 0;
cposition = get(handles.hera_redprogressbar,'Position');

%Set the progress bar text
set(handles.hera_progressbartext,'String','Tracing heartbeat on CWT...');


%Determining heart rate frequency by tracing top of transform for every time-bin
topIndices = nan(1,length(handles.currentfile.matfile.sampleTime)); %Predefine length
for index = startM:endM; %For all time indices
    %index
    if index <= startM+100 %Find starting point
        topIndices(index) = find( handles.currentfile.matfile.CWT(:,index) == max(handles.currentfile.matfile.CWT(:,index)) ); %Find index of maximum value
    else %Consecutive maxima must be found within a constricted frequency range
        try %Try this first
            centreIndex = round( mean( topIndices(index - 99 : index - 1) )); %99 Set middle point
            searchDom = (centreIndex - handles.currentfile.matfile.settings.ballmobilityrange/ ...
                handles.currentfile.matfile.aBin) : (centreIndex + ...
                handles.currentfile.matfile.settings.ballmobilityrange/handles.currentfile.matfile.aBin); %Define domain to be searched for max
            relIndex = find( handles.currentfile.matfile.CWT(searchDom,index) == max(handles.currentfile.matfile.CWT(searchDom,index)) ); %Find index of maximum value within searchDom vector
            topIndices(index) = centreIndex - ((length(searchDom) + 1) / 2) + relIndex; %Index transformation
        catch %If there's an error, do this
            topIndices(index) = find( handles.currentfile.matfile.CWT(:,index) == max(handles.currentfile.matfile.CWT(:,index)) ); %Find index of maximum value
        end
    end
    
    %Update progress bar
    cperc = round(((index-startM) / (endM-startM))*100);
    if cperc > pperc
        pperc = cperc;
        cposition(3) = cperc * handles.hera_def.barscalefactor;
        set(handles.hera_redprogressbar,'Position',cposition);
        drawnow
    end

end

%Convert to actual frequencies
handles.currentfile.matfile.topTrace = nan(1,length(handles.currentfile.matfile.sampleTime)); %Predefine length
for index = startM:endM; 
    handles.currentfile.matfile.topTrace(index) = handles.currentfile.matfile.scalFreq(topIndices(index));
end

%Set the saved setting to "no" to avoid losing changes
handles.currentfile.saved = false;

%Draw and exit
handles = drawgui(handles);
guidata(hObject, handles);



%------------------------------------------------------------------
%Calculate Powerspectrum
%------------------------------------------------------------------

% --- Executes on button press in hera_powerspec.
function hera_powerspec_Callback(hObject, eventdata, handles)
% hObject    handle to hera_powerspec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if ~isempty(handles.currentfile.matfile.topTrace)
    
    %Set the progress bar text
    set(handles.hera_progressbartext,'String','Calculating powerspectrum...');
    
    %Set progress halfway (it is not functional yet)
    cperc = 50;
    cposition = get(handles.hera_redprogressbar,'Position');
    cposition(3) = cperc * handles.hera_def.barscalefactor;
    set(handles.hera_redprogressbar,'Position',cposition);
    drawnow;
    
    ftopTrace = handles.currentfile.matfile.topTrace - nanmean_selfmade(handles.currentfile.matfile.topTrace); %Backup

    %Calculate the sample number to begin from and end at using the
    %startCutout and endCutout settings
    startFnd = find(handles.currentfile.matfile.sampleTime > ...
        handles.currentfile.matfile.settings.startCutout);
    startM = startFnd(1);
    endFnd = find(handles.currentfile.matfile.sampleTime > ...
        (handles.currentfile.matfile.sampleTime(end)- ...
        handles.currentfile.matfile.settings.endCutout));
    endM = endFnd(1);

    %Calculate power spectrum
    %!check these settings: may not work for other samplerates!!!
%    [handles.currentfile.matfile.spectrumPower,...
 %       handles.currentfile.matfile.spectrumFrequencies] = ...
  %      pmusic(ftopTrace(startM:endM),500, ...
   %     2^14,handles.currentfile.matfile.settings.samplerate);

    freqVector = 0.05:0.001:0.5;
    [handles.currentfile.matfile.spectrumPower,...
        handles.currentfile.matfile.spectrumFrequencies] = ...
        pwelch(ftopTrace(startM:endM), [], [], freqVector, ...
        handles.currentfile.matfile.settings.samplerate);

    
    
%testing with pre-ibichannel
%    [handles.currentfile.matfile.spectrumPower,...
%        handles.currentfile.matfile.spectrumFrequencies] = ...
%        pmusic(handles.currentfile.matfile.preibichannel(startM:endM),500, ...
%        2^14,handles.currentfile.matfile.settings.samplerate);
   
    handles = CalcOutmeasures(handles);
    handles = drawgui(handles);
    guidata(hObject, handles);

end


%Calculate output measures
function handles = CalcOutmeasures(handles)


%Find indices (in the spectrum) of frequencies at boundaries of selected frequency bands
fFnd = find(handles.currentfile.matfile.spectrumFrequencies > ...
    handles.currentfile.matfile.settings.lowerboundLF);
handles.currentfile.matfile.lowerboundLFIndex = fFnd(1);

fFnd = find(handles.currentfile.matfile.spectrumFrequencies > ...
    handles.currentfile.matfile.settings.upperboundLF);
handles.currentfile.matfile.upperboundLFIndex = fFnd(1);

fFnd = find(handles.currentfile.matfile.spectrumFrequencies > ...
    handles.currentfile.matfile.settings.lowerboundHF);
handles.currentfile.matfile.lowerboundHFIndex = fFnd(1);

fFnd = find(handles.currentfile.matfile.spectrumFrequencies > ...
    handles.currentfile.matfile.settings.upperboundHF);
handles.currentfile.matfile.upperboundHFIndex = fFnd(1);

    
%Calculate power in selected frequency bands
handles.currentfile.matfile.outmeasures.LFPower = sum(handles.currentfile.matfile.spectrumPower(...
    handles.currentfile.matfile.lowerboundLFIndex:handles.currentfile.matfile.upperboundLFIndex));
handles.currentfile.matfile.outmeasures.HFPower = sum(handles.currentfile.matfile.spectrumPower(...
    handles.currentfile.matfile.lowerboundHFIndex:handles.currentfile.matfile.upperboundHFIndex));
handles.currentfile.matfile.outmeasures.LFHFratio = handles.currentfile.matfile.outmeasures.LFPower / ...
    handles.currentfile.matfile.outmeasures.HFPower;
handles.currentfile.matfile.outmeasures.postBPM = nanmean_selfmade(handles.currentfile.matfile.topTrace)*60;



function hera_preBPM_Callback(hObject, eventdata, handles)
% hObject    handle to hera_preBPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_preBPM as text
%        str2double(get(hObject,'String')) returns contents of hera_preBPM as a double


% --- Executes during object creation, after setting all properties.
function hera_preBPM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_preBPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hera_prerMSSD_Callback(hObject, eventdata, handles)
% hObject    handle to hera_prerMSSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_prerMSSD as text
%        str2double(get(hObject,'String')) returns contents of hera_prerMSSD as a double


% --- Executes during object creation, after setting all properties.
function hera_prerMSSD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_prerMSSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hera_postBPM_Callback(hObject, eventdata, handles)
% hObject    handle to hera_postBPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_postBPM as text
%        str2double(get(hObject,'String')) returns contents of hera_postBPM as a double


% --- Executes during object creation, after setting all properties.
function hera_postBPM_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_postBPM (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hera_LFHFrat_Callback(hObject, eventdata, handles)
% hObject    handle to hera_LFHFrat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_LFHFrat as text
%        str2double(get(hObject,'String')) returns contents of hera_LFHFrat as a double


% --- Executes during object creation, after setting all properties.
function hera_LFHFrat_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_LFHFrat (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hera_LFpower_Callback(hObject, eventdata, handles)
% hObject    handle to hera_LFpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_LFpower as text
%        str2double(get(hObject,'String')) returns contents of hera_LFpower as a double


% --- Executes during object creation, after setting all properties.
function hera_LFpower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_LFpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function hera_HFpower_Callback(hObject, eventdata, handles)
% hObject    handle to hera_HFpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_HFpower as text
%        str2double(get(hObject,'String')) returns contents of hera_HFpower as a double


% --- Executes during object creation, after setting all properties.
function hera_HFpower_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_HFpower (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in hera_prereject.
function hera_prereject_Callback(hObject, eventdata, handles)
% hObject    handle to hera_prereject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


    
%----------------------------------------------------------------
%Reject a time period for pre-wavelet calculation
%----------------------------------------------------------------

%Should only work when a file is open and the file actually contains
%pre-detected peaks (in a pre-ibichannel)
if handles.currentfile.open &  ~isempty(handles.currentfile.matfile.preibichannel)

    %Check the current zoom factor and give a warning if no zooming was
    %done
    
    if get(handles.axes1,'XLim') == handles.currentfile.matfile.prepulsezoom
 
        %Tell user that the selecting the entire region is not useful
        msgbox('Cannot reject the entire time window. Use the zoom tool to select a period containing bad data and press this button again.');
        
    else
        
        %Reject selected
        
        nRej = length(handles.currentfile.matfile.prereject);
        cRej = get(handles.axes1,'XLim');
        
        if cRej(2) > handles.currentfile.matfile.sampleTime(end)
            cRej(2) = handles.currentfile.matfile.sampleTime(end);
        end
        
        if cRej(1) < handles.currentfile.matfile.sampleTime(1)
            cRej(1) = handles.currentfile.matfile.sampleTime(1);
        end
        
        handles.currentfile.matfile.prereject{nRej+1} = cRej;

        %Add the new rejection rectangle in the figure
        cYLim=get(handles.plothandle_peaks(1),'YData'); %get from peaks
        hold on
        
        %20180201: added transparency to the red rectangles that indicate
        %the rejection period.
        handles.plothandle_rejected(nRej+1) = ...
           patch(...
            [cRej(1) ,cRej(2) ,cRej(2) , cRej(1)],...
            [cYLim(1),cYLim(1),cYLim(2),cYLim(2)],...
            [1,.7,.7]);
        set(handles.plothandle_rejected(nRej+1),'FaceAlpha',.5);
        set(handles.plothandle_rejected(nRej+1),'LineStyle','none');
        hold off
            
%         set(handles.plothandle_rejected(nRej+1),'LineStyle','none');        
%         handles.plothandle_rejected(nRej+1) = ...
%             rectangle('Position',...
%                 [cRej(1),cYLim(1),cRej(2)-cRej(1),cYLim(2)-cYLim(1)]);
%         
%         set(handles.plothandle_rejected(nRej+1),'FaceColor',[1,.7,.7]);
%         set(handles.plothandle_rejected(nRej+1),'LineStyle','none');
%         hold off
        
        %Recalculate all pre-wavelet output measures
        handles = CalcPreOutMeasures(handles);
        set(handles.hera_preBPM, 'String',handles.currentfile.matfile.outmeasures.preBPM);
        set(handles.hera_prerMSSD, 'String',handles.currentfile.matfile.outmeasures.rMSSD);
        set(handles.hera_prertMSSD, 'String',handles.currentfile.matfile.outmeasures.rtMSSD);

        %handles = drawgui(handles);
        guidata(hObject, handles);
        
    end    
end



%----------------------------------------------------------------
%De-Reject a time period for pre-wavelet calculation
%----------------------------------------------------------------


% --- Executes on button press in hera_clearprereject.
function hera_clearprereject_Callback(hObject, eventdata, handles)
% hObject    handle to hera_clearprereject (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


%Should only work when a file is open and the file actually contains
%pre-detected peaks (in a pre-ibichannel)
if handles.currentfile.open &  ~isempty(handles.currentfile.matfile.preibichannel)

    %Check if aything has been rejected. If not, give a warning.
    nRej = length(handles.currentfile.matfile.prereject);
    
    if nRej ==0
 
        %Tell user that nothing was rejected
        msgbox('Nothing was rejected, so there is nothing to remove.');
        
    else
        
        %Get the rejection closest to the current position
        %it has to be in the figure
        cZoom = get(handles.axes1,'XLim');
        
        %Loop over rejections and determine if they're on screen
        for crej=1:numel(handles.currentfile.matfile.prereject)
            crejtime=handles.currentfile.matfile.prereject{crej};
            
            crejinzoom(crej)= ~(crejtime(2)<cZoom(1) || crejtime(1)>cZoom(2)); %outside window
            crejmid(crej)=mean(crejtime);
            crejdistmid(crej)=abs(mean(crejtime)-mean(cZoom));
        end
        
        %get the mininum distance one, check if it's on screen, then go
        [irr,cderej]=min(crejdistmid);
        
        if crejinzoom(crejinzoom)~=1
            msgbox('Unclear what to rejection to remove.');
        end
                
        %De-reject the chosen rejection
        newRej = ones(1,numel(handles.currentfile.matfile.prereject));
        newRej(cderej)=0;
        handles.currentfile.matfile.prereject= ...
            handles.currentfile.matfile.prereject(find(newRej==1));
        
        %Also remove from the figure
        delete(handles.plothandle_rejected(cderej));
        handles.plothandle_rejected = ...
            handles.plothandle_rejected(find(newRej==1));
        
        %Redo the calculations
        handles = CalcPreOutMeasures(handles);
        set(handles.hera_preBPM, 'String',handles.currentfile.matfile.outmeasures.preBPM);
        set(handles.hera_prerMSSD, 'String',handles.currentfile.matfile.outmeasures.rMSSD);
        set(handles.hera_prertMSSD, 'String',handles.currentfile.matfile.outmeasures.rtMSSD);
        
        %Finish
        guidata(hObject, handles);
        
    end    
end


%----------------------------------------------------------------
%Calculate pre-wavelet output measures
%----------------------------------------------------------------

function handles = CalcPreOutMeasures(handles)

%Reject all beats within rejected areas by looping over rejected periods
ibis = handles.currentfile.matfile.preibitimeseries;
rejibis = false(1,length(ibis));
ibistarttimes = handles.currentfile.matfile.prepeakTimes(1:end-1);
ibiendtimes = handles.currentfile.matfile.prepeakTimes(2:end);

for rej = 1:length(handles.currentfile.matfile.prereject)
    
    crej = handles.currentfile.matfile.prereject{rej};
    cibistart = find(ibistarttimes > crej(1));
    cibistart = cibistart(1);
    
    cibiend = find(ibiendtimes > crej(2));
    if isempty(cibiend)
        cibiend = numel(ibiendtimes)+1;
    end
    cibiend = cibiend(1)-1;
    rejibis(cibistart:cibiend) = true;
    
end

%Calculate resulting BPM after rejections
%and update settings for tracing after wavelet
handles.currentfile.matfile.settings.avgHeartbeat = 1000/mean(ibis(~rejibis));
%and also set the output measure pre-BPM
%handles.currentfile.matfile.outmeasures.preBPM = 60000/mean(ibis(~rejibis));
handles.currentfile.matfile.outmeasures.preBPM = numel(ibis(~rejibis)) / (sum(ibis(~rejibis))/60000);



%Now loop over all non-rejected periods to calculate rMSSD
rejibis = [true,rejibis,true];
periodstarts = find(diff(rejibis) == -1);
periodends = find(diff(rejibis) == 1)-1;

%loop over periods and fill up successive differences vector
SDs = [];
for period =1:length(periodstarts)
    cperiodibis = ibis(periodstarts(period):periodends(period));
    
    %add this period to SDs vector;
    SDs=[SDs,diff(cperiodibis)];
    
end
   
%Now Calculate rMSSD from SDs
handles.currentfile.matfile.outmeasures.rMSSD = sqrt(mean(SDs.^2));

%20200602: Now also calculate rtMSSD from SDs, using truncated means
%This uses trimmean with 20% of data removed, ie it removes the upper and
%lower 10% of the distribution.
SSDs=SDs.^2;
tMSSDs=trimmean(SSDs,20);
rtMSSD=sqrt(tMSSDs);
handles.currentfile.matfile.outmeasures.rtMSSD = rtMSSD;


%----------------------------------------------------------------
%Insert a peak in pre-wavelet data
%----------------------------------------------------------------

% --- Executes on button press in hera_insertpeak.
function hera_insertpeak_Callback(hObject, eventdata, handles)
% hObject    handle to hera_insertpeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.currentfile.open &  ~isempty(handles.currentfile.matfile.preibichannel)
    
    %Check if user has zoomed on a single peak, otherwise warn
    zoomedperiod = get(handles.axes1,'XLim');
    zoomloc = find(handles.currentfile.matfile.sampleTime > zoomedperiod(1) & ...
        handles.currentfile.matfile.sampleTime < zoomedperiod(2));
    zoomeddata = handles.currentfile.matfile.bpf_pulsedata(zoomloc);
    warning off
    zoomedpeaks = findpeaks(zoomeddata);
    warning on
    
    if numel(zoomedpeaks) == 0
        
        answ = questdlg('No peaks found in zoomed selection. Take the middle?');
        if ~strmatch(answ,'Yes')
            return
        end
        newpeakloc = round((zoomloc(end)-zoomloc(1))/2) + zoomloc(1);
        
    elseif numel(zoomedpeaks) > 1
        
        %Behavior changed 20130316. Instead of error message, just take the
        %peak that is closest to the middle of the figure.
        [zoomedpeaks,zoomedpeaklocs] = findpeaks(zoomeddata);
        
        [irr,newpeakloc]=min(abs(zoomedpeaklocs-size(zoomeddata,2)/2));
        newpeakloc=newpeakloc(1);
        newpeakloc=zoomedpeaklocs(newpeakloc);
        newpeakloc=newpeakloc+zoomloc(1)-1;
        
        
    else
        
        newpeakloc = find(zoomeddata==zoomedpeaks)+zoomloc(1)-1;
        
    end
    newpeakTime = handles.currentfile.matfile.sampleTime(newpeakloc);
    newpeakIBIloc = find(handles.currentfile.matfile.prepeakTimes>newpeakTime);
    newpeakIBIloc = newpeakIBIloc(1);
    
    %Give a warning message if there's already a peak there
    if numel(find(handles.currentfile.matfile.prepeaklocs==newpeakloc))>0
        msgbox('There is already a peak in this location.');
        return
    end
    
    %Now put in the peak, and recalculate everything
    handles.currentfile.matfile.prepeakTimes = ...
        [handles.currentfile.matfile.prepeakTimes(1:newpeakIBIloc-1),...
        handles.currentfile.matfile.sampleTime(newpeakloc), ...
        handles.currentfile.matfile.prepeakTimes(newpeakIBIloc:end)];
    handles.currentfile.matfile.prepeaklocs = ...
        [handles.currentfile.matfile.prepeaklocs(1:newpeakIBIloc-1),newpeakloc, ...
            handles.currentfile.matfile.prepeaklocs(newpeakIBIloc:end)];
    
    %Added 20130316: keep the zoomed area
    %handles.keepzoom = [get(handles.axes1,'XLim');get(handles.axes1,'YLim')];

    %Add the new peak in the figure, adjust the figure handles (20130316)
    cYLim=get(handles.plothandle_peaks(1),'YData');
    xPos=handles.currentfile.matfile.sampleTime(newpeakloc);
    hold on
    newpeakplothandle = plot([xPos,xPos],cYLim,'Color',[.8,.8,.8],'linewidth',.5);
    
    %And put in the plot handle in the correct location
    handles.plothandle_peaks = [handles.plothandle_peaks(1:newpeakIBIloc-1),...
        newpeakplothandle,handles.plothandle_peaks(newpeakIBIloc:end)];
        
    %Call function for creation of pre-ibichannel and ibitimeseries
    %from prepeakTimes
    handles = CalcIBIs(handles);
    
    %Update the IBI channel
    delete(handles.plothandle_preibichannel);
    handles.plothandle_preibichannel = ...
        plot(handles.xaxis,handles.currentfile.matfile.preibichannel,'r','linewidth',1);
    hold off
    
    %Recalculate all pre-wavelet output measures
    handles = CalcPreOutMeasures(handles);
    set(handles.hera_preBPM, 'String',handles.currentfile.matfile.outmeasures.preBPM);
    set(handles.hera_prerMSSD, 'String',handles.currentfile.matfile.outmeasures.rMSSD);
    set(handles.hera_prertMSSD, 'String',handles.currentfile.matfile.outmeasures.rtMSSD);

    %handles = drawgui(handles);
    guidata(hObject, handles);
end



%----------------------------------------------------------------
%Remove a peak from pre-wavelet data
%----------------------------------------------------------------

% --- Executes on button press in hera_removepeak.
function hera_removepeak_Callback(hObject, eventdata, handles)
% hObject    handle to hera_removepeak (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


if handles.currentfile.open &  ~isempty(handles.currentfile.matfile.preibichannel)
    
    %Get indices of peaks in current zoom window
    peaktimes = handles.currentfile.matfile.prepeakTimes;
    zoomedperiod = get(handles.axes1,'XLim');
    rejpeakloc = find(peaktimes > zoomedperiod(1) & peaktimes < zoomedperiod(2));
    %Warning if there are no peaks
    if numel(rejpeakloc)==0
        msgbox('No peak found!')
        return
    end
    
    %Calculate the time at the middle of the figure
    midTime = (zoomedperiod(2)-zoomedperiod(1))./2+zoomedperiod(1);
    %Calculate distances to mid time
    abs(peaktimes(rejpeakloc) - midTime);
    [irr,selPeak]=min(abs(peaktimes(rejpeakloc) - midTime));
    %determine the to-be-rejected peak
    rejpeakloc=rejpeakloc(selPeak(1));
        
    %Now remove this peak
    newpeaks = true(length(peaktimes),1);
    newpeaks(rejpeakloc)=false;
    
    handles.currentfile.matfile.prepeakTimes = ...
        handles.currentfile.matfile.prepeakTimes(newpeaks);
    
    handles.currentfile.matfile.prepeaklocs = ...
        handles.currentfile.matfile.prepeaklocs(newpeaks);

    %Remove the peak from the figure
    delete(handles.plothandle_peaks(rejpeakloc));
    
    %Remover the peak from the peak plot handles
    handles.plothandle_peaks=handles.plothandle_peaks(newpeaks);
    
    %Call function for creation of pre-ibichannel and ibitimeseries
    %from prepeakTimes
    handles = CalcIBIs(handles);
    
    %Update the IBI channel in the figure
    delete(handles.plothandle_preibichannel);
    hold on
    handles.plothandle_preibichannel = ...
        plot(handles.xaxis,handles.currentfile.matfile.preibichannel,'r','linewidth',1);
    hold off
        
    %Recalculate all pre-wavelet output measures
    handles = CalcPreOutMeasures(handles);
    set(handles.hera_preBPM, 'String',handles.currentfile.matfile.outmeasures.preBPM);
    set(handles.hera_prerMSSD, 'String',handles.currentfile.matfile.outmeasures.rMSSD);
    set(handles.hera_prertMSSD, 'String',handles.currentfile.matfile.outmeasures.rtMSSD);
    
    %handles = drawgui(handles);
    guidata(hObject, handles);
end


%----------------------------------------------------------------
%Calculate IBI timeseries and channel
%----------------------------------------------------------------

function handles = CalcIBIs(handles)


%Loop over all peaks, if any
if ~isempty(handles.currentfile.matfile.prepeakTimes)
    
    %Get pre-peak locations
    prepeaklocs = handles.currentfile.matfile.prepeaklocs;
    
    %Calculate and store the pre-IBI timeseries
    handles.currentfile.matfile.preibitimeseries = ...
        diff(handles.currentfile.matfile.prepeakTimes).*1000;

    for cpeak = 1:length(prepeaklocs)

        %Create the pre-IBIchannel
        if cpeak>1
            handles.currentfile.matfile.preibichannel(prepeaklocs(cpeak-1):prepeaklocs(cpeak)) = ...
                handles.currentfile.matfile.preibitimeseries(cpeak-1);
        end
    end

    %Now fill up the start and end
    if prepeaklocs(1)>1 
        handles.currentfile.matfile.preibichannel(1:prepeaklocs(1)-1) = ...
            handles.currentfile.matfile.preibitimeseries(1);
    end

    if prepeaklocs(length(prepeaklocs))<length(handles.currentfile.matfile.rawpulsedata)
        handles.currentfile.matfile.preibichannel(prepeaklocs(length(prepeaklocs)):length(handles.currentfile.matfile.rawpulsedata)) = ...
            handles.currentfile.matfile.preibitimeseries(length(handles.currentfile.matfile.preibitimeseries));
    end
end


        


% --- Executes on button press in hera_collectoutput.
function hera_collectoutput_Callback(hObject, eventdata, handles)
% hObject    handle to hera_collectoutput (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

%Check if SPM is available, if so, use SPM's file selector.
if exist('spm_select') == 2
    ffiles = spm_select(inf,'any','Select processed .puls files',[],pwd,['^*.',handles.hera_def.pulsefile_extension]);
    %If cancelled
    if numel(ffiles)==0
        return
    end
else
    %Otherwise, use another file selector
    answ = 'Yes';
    nDirs = 1;
    ffiles = [];
    while numel(strmatch(answ,'Yes'))==1
        %Open a dialog box for opening multiple files
        [pulsfilenames{nDirs},pulsdir{nDirs}]=uigetfiles('Select files to collect processed data from.');
        
        for i=1:size(pulsfilenames{nDirs},1)
           ffiles = strvcat(ffiles,[ pulsdir{nDirs},pulsfilenames{nDirs}(i,:)]);
        end
        
        nDirs=nDirs+1;
        drawnow;
        answ = questdlg('Add more files (e.g., from a different directory)?');
        drawnow;
    end
end

%Create top line with variable names
outascii = [double('filepath'),9,double('filename'),9,double('IBI_BPM'),9,double('rMSSD'),9,double('rtMSSD'),9,double('wavBPM'),9,double('LFpower'),9,double('HFpower'),9,double('LFHFratio'),13,10];
%Load all files and collect the output data
for i=1:size(ffiles,1)

    cfile=ffiles(i,:);
    extloc = strfind(cfile,['.',handles.hera_def.pulsefile_extension]);
    cmatfile = [cfile(1:extloc(length(extloc))),'mat'];
    if ~exist(cmatfile)
        msgbox(['No ".mat" file found associated with the file: ',cfile,'. Process this file first before collecting output.']);
        return
    else
        cmatfilecontent = load(cmatfile);
        
        %Update 20200415: make sure collect for output works with missing
        %fields, including rtMSSD
        %For each field, check if it exists and if not make an empty one
        if ~isfield(cmatfilecontent.matfile.outmeasures,'preBPM')
            cmatfilecontent.matfile.outmeasures.preBPM=[];
        end
        if ~isfield(cmatfilecontent.matfile.outmeasures,'rMSSD')
            cmatfilecontent.matfile.outmeasures.rMSSD=[];
        end
        if ~isfield(cmatfilecontent.matfile.outmeasures,'rtMSSD')
            cmatfilecontent.matfile.outmeasures.rtMSSD=[];
        end
        if ~isfield(cmatfilecontent.matfile.outmeasures,'postBPM')
            cmatfilecontent.matfile.outmeasures.postBPM=[];
        end
        if ~isfield(cmatfilecontent.matfile.outmeasures,'LFPower')
            cmatfilecontent.matfile.outmeasures.LFPower=[];
        end
        if ~isfield(cmatfilecontent.matfile.outmeasures,'HFPower')
            cmatfilecontent.matfile.outmeasures.HFPower=[];
        end
        if ~isfield(cmatfilecontent.matfile.outmeasures,'LFHFratio')
            cmatfilecontent.matfile.outmeasures.LFHFratio=[];
        end
        
        [cpath,cfile,cext]=fileparts(cmatfile);
        
        %Copy the data into the file, and leave empty values empty    
        outascii = [outascii, double(cpath), 9,double(cfile),9];
        outascii = [outascii, double(num2str(cmatfilecontent.matfile.outmeasures.preBPM)), 9];
        outascii = [outascii, double(num2str(cmatfilecontent.matfile.outmeasures.rMSSD)), 9];
        outascii = [outascii, double(num2str(cmatfilecontent.matfile.outmeasures.rtMSSD)), 9];
        outascii = [outascii, double(num2str(cmatfilecontent.matfile.outmeasures.postBPM)), 9];
        outascii = [outascii, double(num2str(cmatfilecontent.matfile.outmeasures.LFPower)), 9];
        outascii = [outascii, double(num2str(cmatfilecontent.matfile.outmeasures.HFPower)), 9];
        outascii = [outascii, double(num2str(cmatfilecontent.matfile.outmeasures.LFHFratio)), 13, 10];
      
        clear cmatfilecontent
    end
end

filechosen=false;
while ~filechosen
    [outfilename,outfilepath,answ] = uiputfile('Ascii text file, *.txt');
    
    %If user clicked cancel
    if outfilename ==0 & outfilepath ==0
        return
    end
    
    fileopened=false;
    while ~fileopened
        fid=fopen(fullfile(outfilepath,outfilename),'w');
        if fid ==-1
           ans = questdlg('Cannot open file for writing. Try again?');

           if ~strcmp(ans,'Yes')
               fileopened = true; %ie choose another file because filechosen = false
           end
        else
            fileopened = true;
            filechosen = true;
        end
    end
end

fwrite(fid,outascii);
fclose(fid);




%----------------------------------------------------------------------
function [filenames, pathname] = uigetfiles(DialogTitle)
% This is a Java interfaced version of UIGETFILES, that brings multiple file
% open dialog box. 
%
% [filenames, pathname] = uigetfiles; displays a dialog box file browser
% from which the user can select multiple files.  The selected files are
% returned to FILENAMES as an arrayed strings. The directory containing
% these files is returned to PATHNAME as a string. 
%
% A successful return occurs only if the files exist.  If the user selects
% a  file that does not exist, an error message is displayed to the command
% line.  If the Cancel button is selected, zero is assigned to FILENAMES
% and current directory is assigned to PATHNAME. 
% 
% This program has an equivalent function to that of a C version of
% "uigetfiles.dll" downloadable from www.mathworks.com under support, file
% exchange (ID: 331). 
%
% It should work for matlab with Java 1.3.1 or newer.
%
% Shanrong Zhang
% Department of Radiology
% University of Texas Southwestern Medical Center
% 02/09/2004
%
% email: shanrong.zhang@utsouthwestern.edu
% mainFrame = com.mathworks.ide.desktop.MLDesktop.getMLDesktop.getMainFrame;
%
% Additions: Title change using an argument (EJH)

filechooser = javax.swing.JFileChooser(pwd);
filechooser.setDialogTitle(DialogTitle);
filechooser.setMultiSelectionEnabled(true);
filechooser.setFileSelectionMode(filechooser.FILES_ONLY);
selectionStatus = filechooser.showOpenDialog(com.mathworks.mwswing.MJFrame); 

if selectionStatus == filechooser.APPROVE_OPTION
    pathname = [char(filechooser.getCurrentDirectory.getPath), ...
                java.io.File.separatorChar];
    selectedfiles = filechooser.getSelectedFiles;
    for k = 1:1:size(selectedfiles)
        filenames(k) = selectedfiles(k).getName;
    end
    filenames = char(filenames);  
else
    pathname = pwd;
    filenames = 0;
end




%BandPassFilter
%digital band pass filter
%
%(C) Erno Hermans 2004
%usage: bpf(data,samplerate,highpasscutoff,lowpasscutoff)
%
%data: any vector
%samplerate in Hertz
%2 cutoffs in Hertz
%
%use cutoff below zero to turn off (e.g., no highpassfilter)

function out=hera_bpf(data,samplerate,hpcutoff,lpcutoff);


%------------------------ Band Pass Filter ----------------------------



duration = length(data)/samplerate; %duration of input in sec
    
for teller=1:(length(data)/2);
   hertzen(teller+1)= 1/(duration/teller);
   hertzen(length(data)+1-teller)= 1/(duration/teller);
end

fftdata=fft(data);              %fourier transformed data

if lpcutoff>0
    fftdata(hertzen>lpcutoff)=0;    %apply low pass filter
end

if hpcutoff>0
    fftdata(hertzen<hpcutoff)=0;    %apply high pass filter
end

out=real(ifft(fftdata));


% --- Executes on button press in hera_exportIBIs.
function hera_exportIBIs_Callback(hObject, eventdata, handles)
% hObject    handle to hera_exportIBIs (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)




%Loop over all peaks, if any
if handles.currentfile.open   
    if ~isempty(handles.currentfile.matfile.prepeakTimes)

        if ~isempty(handles.currentfile.matfile.prereject)
           uiwait(msgbox('Some parts of the signal were rejected. Note that signal is exported including rejected periods.','Warning','modal'));
        end
        
        %Use a save as type dialog
        [outfilename,outfilepath,answ] = uiputfile([fullfile(handles.currentfile.pulsdir,handles.currentfile.matfilename(1:end-4)),'.txt'], ...
            'Save IBI time series data to ascii text file');
        
        %If user clicked cancel
        if outfilename ==0 & outfilepath ==0
            return
        end
        
        %Save the pre-ibi timeseries
        preibitimeseries = handles.currentfile.matfile.preibitimeseries;
        save(fullfile(outfilepath,outfilename),'preibitimeseries','-ascii');
    end
end


%nanmean selfmade, does exactly the same as nanmean, but doesn't require a
%statistics toolbox license
function outdata = nanmean_selfmade(indata)

outdata = mean(indata(~isnan(indata)));



function hera_prertMSSD_Callback(hObject, eventdata, handles)
% hObject    handle to hera_prertMSSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of hera_prertMSSD as text
%        str2double(get(hObject,'String')) returns contents of hera_prertMSSD as a double


% --- Executes during object creation, after setting all properties.
function hera_prertMSSD_CreateFcn(hObject, eventdata, handles)
% hObject    handle to hera_prertMSSD (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
