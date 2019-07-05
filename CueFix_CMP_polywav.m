%% CueFix_FISHmod_polywav
% Description: 
%     For marker fixing using when data has several distinct wav files
%     Loads biosemi data from .bdf file using the biosig toolbox. Compares
%     markers of imported data to embedded markers in .wav file. 
%     Parses the Presentation logfile to get other cues, like attention
%     cues or other stimuli that might have been lost due to port conflict

% INPUTS:
%     projdir: path to project directory
%     datfiles: cell array, columns are filenames involved in makrer fixing for this dataset 
%         col 1: string, .bdf data filename
%         col 2: string, corresponding .wav file used in that dataset
%         col 3: string, logfile name
%     refchan: double, channel indices to be used as a reference in pop_biosig,
%       when more than one channel is specified, will be referenced to avg
%       of thos channels
%     eeg_save: logical, saves a .set file with correct EEG.event and
%       urevent structure
%     yplot: logical, creates plot showing covariance results 

% VERY IMPORTANT PARAMETERS TO CHECK BEFORE RUNNING: 
%     nwavs: number of wav files that 'findpeaks' will look for in covariance function
%     datdir: directory that has .bdf file, and where you will save your new .set file
%     logdir: directory with Presentation logfile
%     line with pop_biosig command, make sure that reference channels and
%     import channels are all correct. 

% Andrew Kessler 05/14/2015
% Miller Lab, CMB, UC Davis

function CueFix_CMP_polywav(projdir,datfile,eeg_save,yplot); 

if nargin<2
    eeg_save = 1; % default save new .set file with fixed events structure
    yplot = 1; % default is to plot
end

datdir = fullfile(projdir,'Raw Data'); dircheck(datdir);
outdir = fullfile(datdir,'Marker Fixed'); dircheck(outdir); 
logdir = fullfile(projdir,'Logs'); dircheck(logdir);
wavdir = fullfile(projdir,'stim'); dircheck(wavdir)

%% Parse logfile

[MrkrList,PosList,UncList,vidnames,vidlines,sndlines,plines,plat] = log_parse(datfile{3},logdir); 

ainds = find(strncmp(MrkrList,'video_',6)); 

nwavs = length(sndlines); 

ainds = find(strncmp(MrkrList,'video_',6)); % find new video event indices, after editing out pauses above... 
acuestmp = MrkrList(ainds); % attn cues
alats = cell2mat(cellfun(@str2num,PosList(ainds),'uni',0)); % remember that logfile output is in 10th of ms!! 
auncs = cell2mat(cellfun(@str2num,UncList(ainds),'uni',0)); 
if any(find(auncs>2))
    warning(['Uncertainties for video cues > 0.2 ms in ' datfile{1} ', please inspect logfile!!!!!!!!!!!!!!!!!!!!'])
else
    display(['Uncertainties for video cues < 0.2 ms in ' datfile{1} ', looks good!!!!!!!!!!!!!!!!!!!!!'])
end





%% VERY IMPORTANT TO SPECIFY REFERENCE CHANNEL(S) HERE, AND WHICH CHANNELS
% SHOULD BE IMPORTED
% channels 35 and 36 are earlobes, 33 and 44 are mastoids, 

EEG = pop_biosig(fullfile(datdir,datfile{1}));
FsEEG = EEG.srate; 
setn = regexprep(datfile{1},'.bdf','');
EEG.setname= setn;

colors = {'k', 'r', 'g', 'b', 'c', 'm'}; 
lcolor = repmat(colors,1,4); 
% add different line options to each color instance
for x = 1:length(colors)
    lcolor{length(colors)+x} = [lcolor{length(colors)+x} '--']; 
    lcolor{2*length(colors)+x} = [lcolor{2*length(colors)+x} '-.']; 
    lcolor{3*length(colors)+x} = [lcolor{3*length(colors)+x} '']; 
end
lcolor = {lcolor{:} lcolor{:}};

%% Remove flat channels, identify by variance between 0.002 and 10 sec 
% occasionally files have large offset in first several samples, so avoid
% this
startsamp = round(EEG.srate*0.002); stopsamp = EEG.srate*10; 
datchunk = EEG.data(:,startsamp:stopsamp);
datchunk = datchunk-repmat(mean(datchunk,2),1,size(datchunk,2)); % subtract mean of this segment
varcheck = var(datchunk,0,2);
rmchan = find(varcheck < mean(varcheck)/100); % channels to remove

%% Alternative method of deleting channels, using pwelch of first 10 sec,
% then checking average power for frequencies between 5 and 100 Hz

% for j = 1:EEG.nbchan        
%     [px(j,:),f] = pwelch([EEG.data(j,1:FsEEG*10)-mean(EEG.data(j,1:FsEEG*10))],FsEEG,FsEEG/2,FsEEG,FsEEG);
% %     hold(htmp,'on');
% %     plot(htmp,f(5:FsEEG/16)',px(j,5:FsEEG/16),lcolor{j});
% end
% 
% pavg = mean(px(:,5:100),2);
% rmstart = find(pavg<10*min(pavg),1,'first'); % first blank EXG channel
% rmchan = rmstart:EEG.nbchan; 



%% part where channels are actually removed, commented out for now... 
% if ~isempty(rmchan)
%     display(['Removing channels ' sprintf('%d,',rmchan) ' in ' datfile{1}])
%     EEG = pop_select(EEG,'channel',setdiff([1:EEG.nbchan],rmchan));
%     EEG = eeg_checkset( EEG );
% else
%     warning(['No flat channels removed in ' datfile{1} ', '  num2str(EEG.nbchan) ...
%         ' channels retained, please verify this is correct']);
% end



% Marker info from EEG data
markers = cell2mat({EEG.urevent(:).type}); 
datlate = cell2mat({EEG.urevent(:).latency}); 


%% Uncomment if you want to find compression ratio using iterative method

% % Marker info from .wav file
% wavnames = datfile{2}{1}; 
% pbmrk = datfile{2}{end};
% carr = []; cuelabarr = {}; cuelatearr = {};
% for m = 1:length(wavnames)
%     [y,FsWAV] = audioread([wavdir wavnames{m}]); 
%     nwavsamps = length(y); clear y; 
%     [cuelat cuelab]=readWavCue([wavdir wavnames{m}]);
% 
%     cuelat = double(cuelat) + 1; % puts first marker at zero sample, so shift all by 1 sample... 
%     cuelab = cuelab;  
%     cuelatwav = double(cuelat);
%     cuelab = cell2mat(cellfun(@str2num,cuelab,'uni',false)); 
% 
%     %% Find compression ratio between EEG data and .wav file data 
%     Rfs = FsEEG/FsWAV;
%     Rcits = 4;
%     Rcdeltastart = [-1e-4:1e-5:1e-4];
%     for g = 1:Rcits
%         clear Rcdelta
%         if g == 1
%             Rcdelta = Rcdeltastart;
%         else
%             Rcmax = Rcdeltaold(maxind);
%             if maxind == 1 || maxind == length(Rcdeltaold)
%                 Rcdelta = Rcmax + Rcdeltastart;
%             else
%                 Rcdelta  = Rcmax + (Rcdeltastart/10);
%             end    
%         end
%         % Rcdelta = [-6e-6:5e-8:-5e-6];
%         maxpks = [];
%         for h = 1:length(Rcdelta)
%             Rc = Rfs + Rcdelta(h);
% %     load([projdir 'Rc.mat'])
% %     if abs(Rfs-Rc)/Rfs*100 > 2
% %         warning('Compression ratio is very different than ratio of sampling rates, covariance probably wont work!!!!!!!!!');
% %     end
%             cuelate = round(cuelat*Rc); 
% 
%             %% Lump covariance, using all cues in 1 pass
% 
%             uncues = unique(cuelab); 
%             unmarks = unique(markers); 
%             unmarks(find(unmarks<10)) = [];
%             markpat = repmat([-1 1],1,length(unmarks)*2); % make longer than we need, any number x2 is even
% 
%             if pbmrk
%                 markpat = repmat([-1 -1 1 1],1,length(unmarks));
%             end
% 
%             temp = zeros(1,cuelate(end)); 
%             dattemp = zeros(1, datlate(end) - datlate(1) + 1); 
% 
%             for j = 1:length(unmarks)
% 
%                 cueinds = find(cuelab == unmarks(j));
%                 mrkinds = find(markers == unmarks(j));
%                 if ~isempty(cueinds)
%                     temp(cuelate(cueinds)) = markpat(j);
%                 end
%                 if ~isempty(mrkinds)
%                     dattemp(datlate(mrkinds)) = markpat(j); 
%                 end
%             end
% 
%             [c lag] = xcov(dattemp,temp); 
%             maxpks(h) = max(c);
%             if g == Rcits
%                 tmpfig = figure; plot(lag/FsEEG,c);
%                 title(gca,['Rcdelta = ' num2str(Rcdelta(h))]); 
%             end
%         end
%         [maxval maxind] = max(maxpks); % indices into Rcdelta
%         Rcdeltaold = Rcdelta; 
%         
%         display(['Finishing iteration ' num2str(g) ', best Rcdelta= ' num2str(Rcdelta(maxind))]);
%         
%     end
%     Rcmax = Rcdeltaold(maxind);
%     Rc = Rfs + Rcmax; 
%     save([projdir 'Rc.mat'],'Rc');
%     % stop debugger here if making compression ratio
%     
%     
%     
%     %% take sum of results
%     carr(m,:) = c; 
%     cuelabarr{m} = cuelab; 
%     cuelatearr{m} = cuelate; 
%     
% end


% Marker info from .wav file
wavnames = datfile{2}{1}; 
pbmrk = datfile{2}{end};
carr = []; cuelabarr = {}; cuelatearr = {};
pks = []; locs = [];
for m = 1:nwavs
    [y,FsWAV] = audioread(fullfile(wavdir, wavnames{m})); 
    nwavsamps = length(y); clear y; 
    [cuelat cuelab]=readWavCue(fullfile(wavdir, wavnames{m}));

    cuelat = double(cuelat) + 1; % puts first marker at zero sample, so shift all by 1 sample... 
    cuelab = cuelab;  
    cuelatwav = double(cuelat);
    cuelab = cell2mat(cellfun(@str2num,cuelab,'uni',false)); 

    %% Find compression ratio between EEG data and .wav file data 
    % Rewrite to use first audio cue in markers, instead of first audio cue in
    % wav cues!!!
    % Maybe make new function, that finds first marker outside of some list. 
    
    % USING FIRST AUDIO CUE
    
    % XXX uncomment line below when not using for DC
%     maybewavons = find(markers<20); % find wav file onsets in data
%     % maybewavons = 1; 
%     wavonsdat = find(markers == mode(markers(maybewavons)));
%     firstmrk = markers(find(markers(wavonsdat(1):end)>90,1,'first')+wavonsdat(1)-1); % find next marker, which is an auditory cue
%     firstcueind = find(cuelab == firstmrk,1,'first'); % find first start marker in wav file
%     firstlat = cuelat(firstcueind); 
%     datstart = datlate(find(markers == firstmrk,1,'first')); 
%     
%     % find marker before next wav file cue, then find that marker in cuelab
%     gapinds = find(diff(datlate) > (round(1.5*FsEEG)));
%     if markers(gapinds(1)) == 20
%         gapinds(1) = gapinds(1)- 1; % pick previous marker
%     end
%     
%     lastmrk = markers(gapinds(1)); 
%     datend = datlate(gapinds(1)); 
%     lastcueind = find(cuelab == lastmrk,1,'last');
%     lastlat = cuelat(lastcueind);
%     
%     display(['Successfully found start of 1st wav file, using marker ' num2str(firstmrk)]); 
%     display(['Successfully found end of 1st wav file, using marker ' num2str(lastmrk)] ); 

    Rfs = FsEEG/FsWAV; 
    load([projdir 'Rc.mat'])
%     Rc = (datend-datstart)/(lastlat-firstlat); 
    if abs(Rfs-Rc)/Rfs*100 > 2
        warning('Compression ratio is very different than ratio of sampling rates, covariance probably wont work!!!!!!!!!');
    end
    cuelate = round(cuelat*Rc); 

    %% Lump covariance, using all cues in 1 pass

    uncues = unique(cuelab); 
    unmarks = unique(markers); 
    unmarks(find(unmarks<10)) = [];
    markpat = repmat([-1 1],1,length(unmarks)*2); % make longer than we need, any number x2 is even

    if pbmrk
        markpat = repmat([-1 -1 1 1],1,length(unmarks));
    end

    temp = zeros(1,cuelate(end)); 
    dattemp = zeros(1, datlate(end) - datlate(1) + 1); 

    for j = 1:length(unmarks)

        cueinds = find(cuelab == unmarks(j));
        mrkinds = find(markers == unmarks(j));
        if ~isempty(cueinds)
            temp(cuelate(cueinds)) = markpat(j);
        end
        if ~isempty(mrkinds)
            dattemp(datlate(mrkinds)) = markpat(j); 
        end
    end

    [c lag] = xcov(dattemp,temp); 
    
    [pks(end+1) locs(end+1)] = max(c);
%     [pks locs] = findpeaks(csum,'NPEAKS',nwavs,'MINPEAKDISTANCE',...
%         round(0.98*nwavsamps/FsWAV*FsEEG),'MINPEAKHEIGHT',0.5*max(csum));

    %% take sum of results
    carr(m,:) = c; 
    cuelabarr{m} = cuelab; 
    cuelatearr{m} = cuelate; 
     
end



csum = sum(carr,1); 
if yplot
    figure; 
    f = plot(lag./FsEEG,csum); hold on; 
    title({'Covariance function for: ' [regexprep(setn,'_',' ') ' using codes: '] num2str(uncues)}); 
    xlabel('lag time (sec)'); 
    xlim([-100 max(lag./FsEEG)]);
    ul = length(uncues);

    if ul > 10
        ncodes = 20;
        nlines = floor(ul/ncodes); % 15 codes per title line
        form = repmat(ncodes,1,nlines); 
        forms = [form (ul-(nlines*ncodes))]; 
        titlemod = mat2cell(uncues,1,forms);
        titlemod = cellfun(@num2str,titlemod,'uni',0); 
        titlemod = {'Covariance function for: ' [regexprep(setn,'_',' ') ' using codes: '] titlemod{:}}; 
        title(titlemod); 
    end
end 

if length(wavnames) ~= nwavs
    warning(['Logfile reports less .wav files present than expected in: ' datfile{1}])
end

wonsets = sort(lag(locs),2,'ascend'); % sample # of wav file onsets
% wonsets(rmwav) = [];
% pks(rmwav) = [];
% locs(rmwav) = [];

if yplot
    plot(lag(locs)./FsEEG,pks,'ro','LineWidth',2)
end


% Assure yourself that these are in fact the onsets, because we are using
% all cues in .wav file for template- so this includes the onset cue.
% This can be confusing because the compression ratio (Rc) was made using
% 1st audio cue (or chirp code etc.)


% get rid of some big stuff for memory constraints
clear dattemp temp c y csum carr

%% Find attention cues from logfile
refmrkrinds = find(strcmp(MrkrList,num2str(cuelabarr{1}(1)))); 

if any([length(refmrkrinds) length(ainds) length(sndlines) length(wonsets)] ~= nwavs)
    warning(['Inconsistent number of videos, sounds, or sound onsets found'...
        'in logfile, please inspect manually ' datfile{3}]);
end


% need to use .wav onset markers in logfile to derive visual onset times
mlist = []; latlist = []; vlats = []; vmrkrs = []; jitset = []; 


for f = 1:length(wonsets)

    audlats = (cuelatearr{f}+wonsets(f))-1; % auditory latencies, shifted by wav onset times from covariance
        
    mlist = [mlist cuelabarr{f}];
    latlist = [latlist audlats];

end



%% Update EEG.event and EEG.urevent structures
eventtmp = struct(); 

% Add a check to make sure none of the events are at latencies beyond file range! 
latcheck = find(latlist>EEG.pnts);
if any(latcheck)
    latlist(latcheck) = []; mlist(latcheck) = []; 
    warning(['Some of new markers are written beyond length of data, removing these markers in ' datfile{1}]);
end

for p = 1:length(mlist)
    eventtmp(p).type = mlist(p); 
    eventtmp(p).latency = latlist(p); 
    eventtmp(p).urevent = p;    
end

EEG.event = eventtmp; 
EEG.urevent = rmfield(eventtmp,'urevent'); 

EEG = pop_editeventvals(EEG,'sort',{'latency' 0}); 
EEG = eeg_checkset(EEG); 

if eeg_save
    % Save a .set file with fixed imported data
    EEG = pop_saveset( EEG, 'filename',[setn '.set'],'filepath',outdir); % ep_ad means epoched and artifact detected
    EEG = eeg_checkset( EEG );
end




end



