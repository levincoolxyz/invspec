function ch = getkeywait(m) 
% GETKEYWAIT - get a key within a time limit
%   CH = GETKEYWAIT(P) waits for a keypress for a maximum of P seconds. P
%   should be a positive number. CH is a double representing the key
%   pressed key as an ascii number, including backspace (8), space (32),
%   enter (13), etc. If a non-ascii key (Ctrl, Alt, etc.) is pressed, CH
%   will be NaN.
%   If no key is pressed within P seconds, -1 is returned, and if something
%   went wrong during excution 0 is returned. 
%
%  See also INPUT,
%           GETKEY (FileExchange)

% author : Jos van der Geest
% email  : jos@jasen.nl

% check input argument
narginchk(1,1);
if numel(m)~=1 || ~isnumeric(m) || ~isfinite(m) || m <= 0,    
    error('Argument should be a single positive number.') ;
end

% set up the timer
tt = timer ;
tt.timerfcn = 'uiresume' ;
tt.startdelay = m ;            

% Set up the figure
% May be the position property should be individually tweaked to avoid visibility
callstr = 'set(gcbf,''Userdata'',double(get(gcbf,''Currentcharacter''))) ; uiresume ' ;
fh = figure(...
    'name','Press a key', ...
    'keypressfcn',callstr, ...
    'windowstyle','modal',... 
    'numbertitle','off', ...
    'position',[0 0  1 1],...
    'userdata',-1) ; 
try
    % Wait for something to happen or the timer to run out
    start(tt) ;    
    uiwait ;
    ch = get(fh,'Userdata') ;
    if isempty(ch), % a non-ascii key was pressed, return a NaN
        ch = NaN ;
    end
catch
    % Something went wrong, return zero.
    ch = 0 ;
end

% clean up the timer ...
stop(tt) ;
delete(tt) ; 
% ... and figure
if ishandle(fh)
    delete(fh) ;
end