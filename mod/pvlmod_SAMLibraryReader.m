function [LibraryDB, ElementNames] = pvlmod_SAMLibraryReader(varargin)
% Imports a library from a text file that matches the following pattern:
%	
% 		library LIBRARYNAME					[1st header line: name]				
% 		type LIBRARYTYPE					[2nd header line: type]
% 		entries 00							[ignored]
% 		entry elementname					[every element entry starts with name]
% 		00:some.lib.name.FIELDX = 00.00		[field and value]
% 		01:some.lib.name.FIELDY = 0.00		[field and value]
% 		...
% 		!									[end of element]
% 		entry anotherelement
% 		00:some.lib.name.FIELDX = 00.00		
% 		01:some.lib.name.FIELDZ = 0.00		[new fields can be introduced at any point]
% 		...
%
% Syntax
%   [LibraryDB, ElementNames] = pvl_SAMLibraryReader()
%   [LibraryDB, ElementNames] = pvl_SAMLibraryReader(LibraryFile)
%
% Description
%
%	Modified from pvl_SAMLibraryReader_CECModules [0]
%
%   pvlmod_SAMLibraryReader reads a System Advisor Model (SAM) library of any type [1]
%   LibraryDB is a vector of structures which describe each 'element' in the SAM library, 
%	ElementNames is a cell column vector of the names of each module in the Library, thus 
%	ElementNames{n} is the same as LibraryDB(n).name. 
%
% Input Parameters:
%   LibraryFile - An optional input string to select which SAM library to
%     read. If omitted, a user dialog box will prompt the user to browse and
%     select a module library. Note that if LibraryFile is input to the
%     function, the standard MATLAB precedence order applies.
%
% Output:
%   LibraryDB - The parameters in the SAM module library. LibraryDB is a column vector of 
%	size [Ne, 1], where each element is a struct with ALL FIELDS LISTED ALONG THE LIBRARY. 
%	And the number of elements Ne is just determined by scaning the library to the end and
%	counting the number of complete entries (finishing with !).
%   ElementNames - A cell array of size [Ne, 1] with the names of the modules in the SAM library. 
%
% Notes:
%    If this function detects that the selected input library is not of a known type, 
%	 (i.e. CECModule, SandiaModule, SandiaInverter, ODM [*]) it will display a warning to the
%	 command window and continue to try and read the library.
%	 [*] Based on the CECModule library, with names changed and (possibly) additional parameters
%		 See CECLibrary_translator().
%
%    Thanks to the PV_LIB team, who wants to thank the SAM team for maintaining the CEC module 
%	 parameter library and allowing for interaction with the library files.
%
% Sources:
%
% [0] PV_LIB v.1.1 web page. http://pvpmc.org
% [1] System Advisor Model web page. https://sam.nrel.gov.

%% Parse the input data
	p = inputParser;
	p.addOptional('LibraryFile', 0, @(x) ischar(x));
	p.parse(varargin{:})

%% Enter Library structure basics
	LibraryHeaderLines = 3; % Number of header lines before beginning data
	NonexistantCharacter = char(172); % A character which is NOT found within the library

	defaultchecker = {'LibraryFile'};
	if any(strcmp(defaultchecker,p.UsingDefaults))
		%% Ask user to get SAM library file
		[FileName, FilePath, FilterIndex] = ...
			uigetfile({'*.samlib;*.modlib','Library Files (*.samlib, *.modlib)';'*.*','All Files'}, ...
				'Select a Library with .samlib or .modlib extension.', 'MultiSelect', 'off');
		if FilterIndex ==0
			error('No .modlib file selected, exiting Library Reader')
		end
		FilePathandName = [FilePath FileName];
	else
		FilePathandName = p.Results.LibraryFile;
	end
	
%% Open the file and read in the header and data separately
	FileID = fopen(FilePathandName);
    assert(FileID > 0,'SAMLibraryReader: Could not open file');
	HeaderDataIn = textscan(FileID, '%s', LibraryHeaderLines, 'Delimiter', NonexistantCharacter);
	RawDataIn = textscan(FileID, '%s', 'Delimiter', NonexistantCharacter);%, 'HeaderLines', LibraryHeaderLines);
	fclose(FileID);
	HeaderDataIn = HeaderDataIn{1};
	RawDataIn = RawDataIn{1};

%% Parse out the library name, library type, and number of entries from the header information
	% This goes with every entry
	LibraryName = textscan(char(HeaderDataIn(1)), '%*8s %s', 'Delimiter', NonexistantCharacter);
	LibraryName = LibraryName{1}{1};

	% This goes with every entry
	LibraryType = textscan(char(HeaderDataIn(2)), '%*s %s', 'Delimiter', ' ', 'MultipleDelimsAsOne', 1);
	LibraryType = LibraryType{1}{1};

	if ~any(strcmp(LibraryType,{'CECModule', 'SandiaModule', 'SandiaInverter', 'ODM','BPD','UF'}))
		warning('pvlmod:SAMLibraryReader:type',...
			'The library type is not recognized. Your data may not be correct.');
	end

	StatedNumberOfEntries = textscan(char(HeaderDataIn(3)), '%*s %f', 'Delimiter',' ', 'MultipleDelimsAsOne', 1);
	StatedNumberOfEntries = StatedNumberOfEntries{1};

%% Scan the library to search for all existing field names, and check the number of entries
	NumberOfEntries = 0;
	fieldnames = cell(40,1); % just to say something
	fieldnames{1} = 'name';
	usednames = 1;
	for j = 1:numel(RawDataIn)
		if strcmp('!',RawDataIn{j})
			NumberOfEntries = NumberOfEntries+1; 
		else
			thisname = textscan(RawDataIn{j},'%02d:%[^.].%[^.].%[^.].%[^=]');
			if isempty(thisname{1}), continue; end
			thisname = thisname{5}{1};
			if ~any(strcmp(thisname,fieldnames))
				usednames = usednames+1;
				fieldnames{usednames} = thisname;
			end
		end
	end
	
	if ~any(strcmp('LibraryName',fieldnames))
		fieldnames{usednames+1} = 'LibraryName';
		usednames = usednames+1;
	end
	if ~any(strcmp('LibraryType',fieldnames))
		fieldnames{usednames+1} = 'LibraryType';
		usednames = usednames+1;
	end
	fieldnames = fieldnames(1:usednames);

	% Create vector of structures with the appropriate fields
	for j = 1:usednames
		emptyStruct.(fieldnames{j}) = [];
	end
	database(1:NumberOfEntries) = emptyStruct;
	database = database(:);
	
	if NumberOfEntries ~= StatedNumberOfEntries
		warning('pvlmod:SAMLibraryReader:Nentries',...
				['The stated number of entries is %d, yet %d were counted.',...
				 'It''s no big deal, but check it.'],StatedNumberOfEntries,NumberOfEntries);
	end

%% Step through every entry and extract the information into the structures
	entry = 1;
	newentry = true;
	for j = 1:numel(RawDataIn)
		if strcmp('!',RawDataIn{j})
			database(entry).LibraryName = LibraryName;
			database(entry).LibraryType = LibraryType;
			entry = entry+1;
			newentry = true;
		else
			if newentry
				thisrow = textscan(RawDataIn{j},'entry %[^\n\r]');
				if ~isempty(thisrow{1})
					database(entry).name = thisrow{1}{1};
					newentry = false;
				else
					continue; % skip it, maybe it's just a blank line
				end
			else
				thisrow = textscan(RawDataIn{j},'%02d:%[^.].%[^.].%[^.].%[^=]=%[^\n\r]');
				if isempty(thisrow{1}), continue; end % skip it
				fieldname = thisrow{5}{1};
				fieldvalue = thisrow{6}{1};
				if isnan(str2double(fieldvalue))
					database(entry).(fieldname) = fieldvalue;
				else
					database(entry).(fieldname) = str2double(fieldvalue);
				end
			end
		end
	end
 
%% Write the information from the CEC Library to the output variables.
	LibraryDB = database;
	ElementNames = cell(numel(LibraryDB,1));
	[ElementNames{:}] = LibraryDB.name;
    
%% MHA: Fix/adjust inverter fields for compatibility with ODM-based models

    if strcmpi(LibraryDB.LibraryType,'SandiaInverter')
        LibraryDB.fPV = @(Pdc,Vdc) pvl_snlinverter(LibraryDB, Vdc, Pdc);
        LibraryDB.TPLim = @(Ta) repmat(LibraryDB.Pac0,size(Ta)); % SANDIA model has no Temp. derate

        % Fix fieldname issues with SAMLIB...
        LibraryDB = replacefield(LibraryDB,'Idcmax','IdcMax');
        LibraryDB.PacMax = LibraryDB.Pac0;
    end
end

function S = replacefield(S,f,g)
    if ~isfield(S,f), return; end
    S.(g) = S.(f);
    S = rmfield(S,f);
end