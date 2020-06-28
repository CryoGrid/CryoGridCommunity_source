%==========================================================================
%  Author: Andriy Nych ( nych.andriy@gmail.com )
% Version:        733341.4155741782200
%==========================================================================
% 
% INI = ini2struct(FileName)
% 
% This function parses INI file FileName and returns it as a structure with
% section names and keys as fields.
% 
% Sections from INI file are returned as fields of INI structure.
% Each fiels (section of INI file) in turn is structure.
% It's fields are variables from the corresponding section of the INI file.
% 
% If INI file contains "oprhan" variables at the beginning, they will be
% added as fields to INI structure.
% 
% Lines starting with ';' and '#' are ignored (comments).
% 
% See example below for more information.
% 
% Usually, INI files allow to put spaces and numbers in section names
% without restrictions as long as section name is between '[' and ']'.
% It makes people crazy to convert them to valid Matlab variables.
% For this purpose Matlab provides GENVARNAME function, which does
%  "Construct a valid MATLAB variable name from a given candidate".
% See 'help genvarname' for more information.
% 
% The INI2STRUCT function uses the GENVARNAME to convert strange INI
% file string into valid Matlab field names.
% 
% [ test.ini ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
%     SectionlessVar1=Oops
%     SectionlessVar2=I did it again ;o)
%     [Application]
%     Title = Cool program
%     LastDir = c:\Far\Far\Away
%     NumberOFSections = 2
%     [1st section]
%     param1 = val1
%     Param 2 = Val 2
%     [Section #2]
%     param1 = val1
%     Param 2 = Val 2
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% The function converts this INI file it to the following structure:
% 
% [ MatLab session (R2006b) ]~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
%  >> INI = ini2struct('test.ini');
%  >> disp(INI)
%         sectionlessvar1: 'Oops'
%         sectionlessvar2: 'I did it again ;o)'
%             application: [1x1 struct]
%             x1stSection: [1x1 struct]
%            section0x232: [1x1 struct]
% 
%  >> disp(INI.application)
%                    title: 'Cool program'
%                  lastdir: 'c:\Far\Far\Away'
%         numberofsections: '2'
% 
%  >> disp(INI.x1stSection)
%         param1: 'val1'
%         param2: 'Val 2'
% 
%  >> disp(INI.section0x232)
%         param1: 'val1'
%         param2: 'Val 2'
% 
% ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
% 
% NOTE.
% WhatToDoWithMyVeryCoolSectionAndVariableNamesInIniFileMyVeryCoolProgramWrites?
% GENVARNAME also does the following:
%   "Any string that exceeds NAMELENGTHMAX is truncated". (doc genvarname)
% Period.
% 
% =========================================================================
function Struct = ini2struct_improved(FileName)

% Parses .ini file
% Returns a structure with section names and keys as fields.
%
% Based on init2struct.m by Andriy Nych
% 2014/02/01

f = fopen(FileName,'r');                    % open file
while ~feof(f)                              % and read until it ends
    s = strtrim(fgetl(f));                  % remove leading/trailing spaces
    if isempty(s) || s(1)==';' || s(1)=='#' % skip empty & comments lines
        continue
    end
    if s(1)=='['                            % section header
        Section = genvarname(strtok(s(2:end), ']'));
        Struct.(Section) = [];              % create field
        continue
    end
    
    [Key,Val] = strtok(s, '=');             % Key = Value ; comment
    Val = strtrim(Val(2:end));              % remove spaces after =
    
    if isempty(Val) || Val(1)==';' || Val(1)=='#' % empty entry
        Val = [];
    elseif Val(1)=='"'                      % double-quoted string
        Val = strtok(Val, '"');
    elseif Val(1)==''''                     % single-quoted string
        Val = strtok(Val, '''');
    else
        Val = strtok(Val, ';');             % remove inline comment
        Val = strtok(Val, '#');             % remove inline comment
        Val = strtrim(Val);                 % remove spaces before comment
        
        [val, status] = str2num(Val);       %#ok<ST2NM>
        if status, Val = val; end           % convert string to number(s)
    end
    
    if ~exist('Section', 'var')             % No section found before
        Struct.(genvarname(Key)) = Val;
    else                                    % Section found before, fill it
        Struct.(Section).(genvarname(Key)) = Val;
    end
end
fclose(f);