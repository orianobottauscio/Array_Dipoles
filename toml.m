classdef toml

methods (Static)

% READ parse TOML data from a file
%
%   READ('file.toml') loads the contents of `file.toml` and parses
%   that data into a MATLAB struct.
%
%   See also FILEREAD, TOML.DECODE

function toml_data = read(filename)
  raw_text = fileread(filename);
  toml_data = toml.decode(raw_text);
end

% DECODE convert TOML to native MATLAB datatypes
%
%   DECODE(toml_str) returns the MATLAB representation of the
%   TOML-formatted data in `toml_str`. If it is invalid TOML, an
%   appropriate exception will be raised.
%
%   See also TOML.READ

function obj_out = decode(toml_str)
%% pre-emptive checking
  % split on newlines
  toml_lines = strsplit(toml_str, {'\n', '\r'});
  % throw out comments
  de_commenter = @(elem) deblank(toml.decomment(elem));
  toml_decommented = cellfun(de_commenter, toml_lines, ...
                             'UniformOutput', false);
  % strip out empty lines
  toml_nonempty = toml_decommented(~cellfun(@isempty, toml_decommented));
  % check for invalid lines
  cellfun(@toml.checkline, toml_nonempty);

%% parsing
  obj_out = struct();
  current_line = 1;
  location_stack = {};

  while current_line <= length(toml_nonempty)
    % recognize a section and store it semantically
    n_brackets = toml.is_section(toml_nonempty{current_line});
    if n_brackets
      section_name = strtrim(toml_nonempty{current_line});
      section_name = section_name(n_brackets+1:end-n_brackets);
      location_stack = toml.parsekey(section_name);
      location_stack = toml.adjust_key_stack(obj_out, location_stack);
      % is it a table or an array of tables?
      if n_brackets == 1
        obj_out = toml.set_nested_field(obj_out, location_stack, struct());
      else
        % if it already exists, append
        try
          existing_val = toml.get_nested_field(obj_out, location_stack);
          location_stack{end + 1} = length(existing_val) + 1; %#ok<*AGROW>
          obj_out = toml.set_nested_field(obj_out, location_stack, struct());
        % if not, pre-populate
        catch
          obj_out = toml.set_nested_field(obj_out, location_stack, {struct()});
          location_stack{end + 1} = 1;
        end
      end
      current_line = current_line + 1;
      continue
    end

    % recognize key-value pairs and add them to the struct
    [key, value] = strtok(toml_nonempty{current_line}, '=');
    key_seq = toml.parsekey(key);
    % ensure we have a complete value
    force = false;
    while true
      value_fix = toml.parsevalue(value, force);
      % loop for possible multiline values
      if isempty(value_fix) && ~iscell(value_fix)
        if current_line < length(toml_nonempty)
          current_line = current_line + 1;
          value = sprintf('%s\n%s', value, toml_nonempty{current_line});
        else
          force = true;
        end
      else
        % convert closed but empty string values to empty char
        if numel(value_fix) == 2 && ...
            (isequal(value_fix, '""') || isequal(value_fix, ''''''))
          value_fix = '';
        end
        break
      end
    end
    obj_out = toml.set_nested_field(obj_out, [location_stack, key_seq], value_fix);
    current_line = current_line + 1;
  end
end



% DECOMMENT remove comment from a single line
%
%   DECOMMENT('# this is a comment') returns ''.
%
%   DECOMMENT('key = "value" # comment') returns 'key = "value"'.
function out = decomment(in)
  all_parts = toml.splitby(in, '#', {'"', ''''});
  out = all_parts{1};
end


% SPLITBY split a string while respecting grouping
%
%   SPLITBY(str, '.', {'"', '[]'}) splits `str` on each period which is
%   not enclosed in square brackets or double quotes.

function parsed = splitby(str, to_split_by, to_respect)
  % make convenient structure out of delimiter arguments
  delims = cellfun(@toml.delim_data, to_respect);
  % walk through the string char by char
  for ch = 1:length(str)
    % check each delimiter (set) for matches
    for delim = 1:length(delims)
      if delims(delim).match_begin(str(ch))
        delims(delim).depth = delims(delim).depth + 1;
        continue
      elseif delims(delim).match_end(str(ch))
        delims(delim).depth = delims(delim).depth - 1;
        continue
      end
    end

    % check for splittability
    if str(ch) == to_split_by
      % check each delimiter (set) for lexical closure
      not_ready = zeros(size(delims));
      for delim = 1:length(delims)
        not_ready(delim) = delims(delim).check_in(delims(delim).depth);
      end

      % if no one complains, mark as splittable
      if ~any(not_ready)
        str(ch) = char(0);
      end
    end
  end

  % actually split the thing
  parsed = strsplit(str, char(0));
end

function dd = delim_data(delim)
  % there is always at least one delimiter
  dd.match_begin = @(c) c == delim(1);
  % start at zero immersion
  dd.depth = 0;
  switch length(delim)
    case 1
      % no intuitive meaning to 'end'
      dd.match_end = @(c) false;
      % as long as there are an even number, we're good
      dd.check_in = @(x) mod(x, 2);
    case 2
      % simple enough
      dd.match_end = @(c) c == delim(2);
      % we want no levels of immersion
      dd.check_in = @(x) min(1, x);
  end
end

% CHECKLINE check a line of TOML data for validity
%
%   CHECKLINE(toml_line) tests that line for invalid forms, and raises
%   an exception if one of them is applicable.

function checkline(in)
  % check for unspecified value
  if in(end) == '='
    error('toml:UnspecifiedValue', ...
          'TOML keys must have a corresponding value.')
  end

  % check for empty (bare) key
  if in(1) == '='
    error('toml:EmptyBareKey', ...
          'TOML bare keys must not be empty.')
  end
end


function n_brackets = is_section(str)
  section_regexp = '^\s*\[{1,2}(.+?\.?)+\]{1,2}$';
  section_name = regexp(str, section_regexp, 'ONCE');

  if ~isempty(section_name)
    if str(2) == '['
      n_brackets = 2;
    else
      n_brackets = 1;
    end
  else
    n_brackets = 0;
  end

end

% PARSEKEY convert a string into a struct pointer
%
%   PARSEKEY(str) splits a string semantically on dots (while respecting
%   quotes), then converts each segment of the resulting pointer into a
%   valid field name with a predictable structure.
%
%   See also PARSEVALUE, SPLITBY

function key = parsekey(str)
  % split on dots, if not inside quotes
  key_seq = toml.splitby(str, '.', {'''', '"'});

  % utility for the following
  uo = {'UniformOutput', false};

  % remove quotes
  dequote = @(elem) regexprep(elem, '["'']+', '');
  key_unquoted = cellfun(dequote, key_seq, uo{:});

  % trim leading and trailing space
  key_trimmed = cellfun(@strtrim, key_unquoted, uo{:});

  % sub underscores for spaces
  despace = @(elem) strrep(elem, ' ', '_');
  key_despaced = cellfun(despace, key_trimmed, uo{:});

  % make it a valid name
  fixname = @(elem) matlab.lang.makeValidName(elem, 'prefix', 'f');
  key = cellfun(fixname, key_despaced, uo{:});

end


% PARSEVALUE parse the corresponding MATLAB object out of a TOML value
%
%   PARSEVALUE('') returns an empty string.
%   PARSEVALUE('0b10') returns 2.
%   PARSEVALUE('0o10') returns 8.
%   PARSEVALUE('0x10') returns 16.
%   PARSEVALUE('10') returns 10.
%   PARSEVALUE('"foo"') returns 'foo'.
%   PARSEVALUE('true') returns logical 1.
%   PARSEVALUE('"\n"') returns a newline character.
%   PARSEVALUE('20180625T07:00Z') returns a datetime object with the
%   value June 25, 2018, 7:00AM UTC.
%   PARSEVALUE('["a", "b"]') returns {'a', 'b'}.
%   PARSEVALUE('[1, 2, 3]') returns [1, 2, 3].
%   PARSEVALUE('{key = "value"}') returns struct('key', 'value').
%
%   See also PARSEKEY

function val = parsevalue(str, force)
%% check for noncompletion
  if isempty(str)
    val = '';
    return
  end

%% default fixes
  % default behavior is direct passthrough
  val = str;

  % remove leading equals sign
  if val(1) == '='
    val = val(2:end);
  end

  % trim space on each side
  trimmed_val = strtrim(val);

%%% Check for string, string literal, array, or table before other datatypes
%% strings
  % basic strings
  if trimmed_val(1) == '"'
    % is it multiline and complete?
    if numel(trimmed_val) > 2 && ...
       isequal(trimmed_val(1:3), '"""') && ...
       numel(trimmed_val) > 3 && ...
       isequal(trimmed_val(end-2:end), '"""')
      % remove quotes
      val = trimmed_val(4:end-3);
      % remove leading newline
      if val(1) == newline
        val = val(2:end);
      end
      % trim whitespace for backslashes
      val = regexprep(val, '\\\n\s+', '');

    % is it complete but not multiline?
    elseif trimmed_val(2) ~= '"' && trimmed_val(end) == '"'
      % remove quotes
      val = trimmed_val(2:end-1);

    % newline in string, tell caller the value is incomplete
    elseif force
      error('toml:IncompleteString', ...
            'String without closing quote: %s', trimmed_val)
    else
      % is it complete but empty?
      if isequal(trimmed_val, '""')
        % set to self so caller doesn't find it empty
        val = '""';
      else
        val = '';
      end
      return
    end

    % common post-processing
    % catch invalid escapes
    invalid_esc = regexp(val, '(?<!\\)\\(\\\\)*([^btnfr"\\uU])', 'match');
    if ~isempty(invalid_esc)
      error('toml:InvalidEscapeSequence', ...
            ['Invalid escape sequence: \', invalid_esc{:}, ...
            '\nFound in this value: ', strrep(str, '\', '\\')])
    end
    % unicode points (only 4 digit hex codes will work in MATLAB)
    ucode_match = '(?<!\\)\\(u[A-Fa-f0-9]{1,4}|U[A-Fa-f0-9]{1,8})';
    ucode_replace = '${char(hex2dec($1(2:end)))}';
    val = regexprep(val, ucode_match, ucode_replace);
    % escaped characters
    val = regexprep(val, '(\\[btnfr"\\])', '${sprintf($1)}');
    return
  end

  % literal strings
  if trimmed_val(1) == ''''
    % is it multiline and complete?
    if numel(trimmed_val) > 2 && ...
       isequal(trimmed_val(1:3), '''''''') && ...
       numel(trimmed_val) > 3 && ...
       isequal(trimmed_val(end-2:end), '''''''')
      % remove quotes
      val = trimmed_val(4:end-3);
      % remove leading newline
      if val(1) == newline
        val = val(2:end);
      end

    % is it complete but not multiline?
    elseif trimmed_val(2) ~= '''' && trimmed_val(end) == ''''
      % remove quotes
      val = trimmed_val(2:end-1);

    % newline in string, tell caller the value is incomplete
    elseif force
      error('toml:IncompleteString', ...
        'String without closing quote: %s', trimmed_val)
    else
      % is it complete but empty?
      if isequal(trimmed_val, '''''')
        % set to self so caller doesn't find it empty
        val = '''''';
      else
        val = '';
      end
    end

    return
  end

%% arrays
  % is it an array?
  if trimmed_val(1) == '['
    % get starting and ending brackets
    beginning_brackets = regexp(trimmed_val, '^\s*\[+[^0-9a-zA-Z"-]*', 'match');
    if isempty(beginning_brackets), beginning_brackets = ''; end
    beginning_brackets = strjoin(split(beginning_brackets), '');
    closing_brackets = regexp(trimmed_val, '[^0-9a-zA-Z"]*\s*\]+$', 'match');
    if isempty(closing_brackets), closing_brackets = ''; end
    closing_brackets = strjoin(split(closing_brackets), ''); %#ok<NASGU>

    % get all opening and closing brackets
    num_opening_brackets = sum(ismember(strjoin(...
      regexp(trimmed_val, '\s*\[+[^0-9a-zA-Z"]*', 'match'), ''), '['));
    num_closing_brackets = sum(ismember(strjoin(...
      regexp(trimmed_val, '\s*\]+[^0-9a-zA-Z"]*', 'match'), ''), ']'));
    
    % is it all here yet?
    if isequal(num_opening_brackets, num_closing_brackets)
      % remove outer brackets
      val = trimmed_val(2:end-1);
      max_dimension = max(size(beginning_brackets));
    else
      if force
        error('toml:IncompleteArray', ...
          'Array without closing bracket: %s', trimmed_val)
      end
      val = [];
      return
    end

    if isempty(val)
      val = {};
      return
    end

    % split array while respecting nesting
    val = toml.splitby(val, ',', {'{}', '[]', '"', ''''});
    val = cellfun(@strtrim, val, 'uniformoutput', false);
    val = val(~cellfun(@isempty, val));
    row_count = sum(cellfun(@(x) isequal(x(1), '[') && ...
      isequal(x(end), ']'), val));
    if row_count == 0, row_count = 1; end
    val = cellfun(@toml.parsevalue, val, 'uniformoutput', false);

    % check homogeneity
    contained_types = cellfun(@class, val, 'uniformoutput', false);
    contained_sizes = cellfun(@numel, val);
    contained_types(strcmp(contained_types, 'double') & contained_sizes > 1) ...
        = deal({'cell'});
    if numel(unique(contained_types)) > 1
      error('toml:HeterogeneousArray', ...
            'All elements of a TOML array must be the same type.')
    elseif all(cellfun(@isnumeric, val))
      val = reshape(val, row_count, []);
      % check if numeric cells have the same number of columns per row
      cell_sizes = cell2mat(cellfun(@size, val, 'UniformOutput', false));
      if numel(unique(cell_sizes(:, 1))) == 1 && numel(unique(cell_sizes(:, 2))) == 1
        % apply dimensions to fully numeric array
        if max_dimension == 2
          max_dimension = 1;
        elseif max_dimension == 1
          max_dimension = 2;
        end
        val = cat(max_dimension, val{:});
      end
    end

    return
  end

%% tables

  if ~isempty(regexp(trimmed_val, '^{(\s*[^=]+\s*=|})', 'ONCE'))
    % is it all here yet?
    if trimmed_val(end) ~= '}'
      if force
        error('toml:IncompleteInlineTable', ...
              'Inline table without closing curly brace: %s', trimmed_val)
      end
      val = [];
      return
    % remove outer brackets
    else
      val = trimmed_val(2:end-1);
    end

    % empty table
    if isempty(val)
      val = struct();
      return
    end

    % split table while respecting nesting
    val = toml.splitby(val, ',', {'{}', '[]', '"', ''''});

    vals = cellfun(@(elem) toml.splitby(elem, '=', {'{}', '[]', '"', ''''}), ...
                   val, 'uniformoutput', false);
    key_names = cellfun(@(elem) toml.parsekey(elem{1}), vals, 'uniformoutput', false);

    val = struct();
    for elem = 1:length(vals)
      val = setfield(val, key_names{elem}{:}, toml.parsevalue(vals{elem}{2}));
    end

    return
  end

%% datetimes

  % make regexes
  is_match = @(s, p) ~isempty(regexp(s, p, 'ONCE'));
  date_regexp = '\d{4}-\d{2}-\d{2}';
  upto24 = ['(', strjoin( ...
      arrayfun(@(elem) sprintf('%02d', elem), 0:23, 'uniformoutput', false), ...
      '|'), ')'];
  fract_sec = '\.\d{1,9}';
  time_regexp = [upto24, ':[0-6]\d:[0-6]\d(', fract_sec, ')?'];
  offset_regexp = ['(Z|[-+]', upto24, ':[0-6]\d(:[0-6]\d)?)$'];

  % see what fits
  has_date = is_match(trimmed_val, date_regexp);
  has_time = is_match(trimmed_val, time_regexp);
  is_datetime = has_date && has_time;
  is_datetime_t = is_datetime && ...
      is_match(trimmed_val, [date_regexp, 'T', time_regexp]);
  has_fr_sec = has_time && is_match(trimmed_val, fract_sec);
  has_offset = has_time && is_match(trimmed_val, offset_regexp);

  % make formats
  date_fmt = 'yyyy-MM-dd';
  time_fmt = 'HH:mm:ss';
  fract_sec_fmt = '.SSSSSSSSS';

  % do what we can with it
  if is_datetime
    dtargs = {};
    if is_datetime_t
      fmt_str = [date_fmt, '''T''', time_fmt];
    else
      fmt_str = [date_fmt, ' ', time_fmt];
    end
    if has_fr_sec
      fmt_str = [fmt_str, fract_sec_fmt];
    end
    if has_offset
      fmt_str = [fmt_str, 'Z'];
      dtargs = {'TimeZone', 'UTC'};
    end
    val = datetime(trimmed_val, 'InputFormat', fmt_str, dtargs{:});
    return
  elseif has_date
    fmt_str = date_fmt;
    val = datetime(trimmed_val, 'InputFormat', fmt_str);
    return
  elseif has_time
    fmt_str = time_fmt;
    if has_fr_sec
      fmt_str = [fmt_str, fract_sec_fmt];
    end
    val = datetime(trimmed_val, 'InputFormat', fmt_str);
    return
  end

%% check for numeric types

  % utils for integers
  is_int = @(t, s, c) all(ismember(t, [s, '_'])) && ...
           length(t) > 3 && isequal(t(1:2), ['0', c]);
  descore = @(t) strrep(t(3:end), '_', '');
  specs.bin = 'b01';
  specs.oct = 'o01234567';
  specs.hex = 'x0123456789abcdefABCDEF';

  dec_num = '([0-9][_0-9]*)?[0-9]';
  dec_int = ['[+-]?', dec_num];
  specs.dec = [ ...
      '^', dec_int ...            % integer part
      '(\.', dec_num, ')?' ...    % fractional part
      '([eE]', dec_int, ')?$' ... % exponential part
              ];

  % binary
  if is_int(trimmed_val, specs.bin, 'b')
    val = bin2dec(descore(trimmed_val));
    return
  % octal
  elseif is_int(trimmed_val, specs.oct, 'o')
    val = base2dec(descore(trimmed_val), 8);
    return
  % hexadecimal
  elseif is_int(trimmed_val, specs.hex, 'x')
    val = hex2dec(descore(trimmed_val));
    return
  % decimal (including floats)
  elseif ~isempty(regexp(trimmed_val, specs.dec, 'ONCE'))
    val = str2double(strrep(val, '_', ''));
    % error for using leading zeros on a decimal integer
    if isfinite(val) && ~mod(val, 1) && val && trimmed_val(1) == '0'
      error('toml:DecIntLeadingZeros', ...
            'Decimal integers may not have leading zeros.')
    end
    return
  end

  % special values of float
  spec_flt = {'inf'; 'nan'};
  op_chars = {'', '+', '-'};
  spec_cmp = reshape( ...
      strcat(repmat(op_chars, 2, 1), repmat(spec_flt, 1, 3)), ...
                     [], 1);
  if any(strcmp(trimmed_val, spec_cmp))
    val = str2double(val);
    return
  elseif any(strcmpi(trimmed_val, spec_cmp))
    error('toml:UppercaseSpecialFloat', ...
          'Special floating-point values must be lowercase.')
  end

%% booleans
  if any(strcmp(trimmed_val, {'true', 'false'}))
    val = ~mod(numel(trimmed_val), 2);
    return
  elseif any(strcmpi(trimmed_val, {'true', 'false'}))
    error('toml:UppercaseBoolean', ...
          'Boolean values must be lowercase.')
  end

%% invalid datatype?
  error('toml:InvalidType', ...
        'Unknown datatype: %s', trimmed_val)
end

% SET_NESTED_FIELD set a value somewhere in a struct
%
%   SET_NESTED_FIELD(obj, indx, val) sets the location denoted by `indx`
%   (a pointer sequence into `obj`) in `obj` equal to `val`, and returns
%   a modified copy of `obj`.
%
%   See also GET_NESTED_FIELD

function obj = set_nested_field(obj, indx, val)
  if length(indx) == 1
    if isstruct(obj)
      if isfield(obj, indx{:})
        switch class(obj.(indx{:}))
          case 'struct'
            if ~isstruct(val) || isempty(fieldnames(val))
              error('toml:RedefinedTable', ...
                    'Tables cannot be redefined.')
            end
          case 'cell'
            if iscell(val) && ...
              ( ...
                isempty(val) || ...
                ( ...
                  isempty(val{1}) || ( ...
                    isstruct(val{1}) && isempty(fieldnames(val{1})) ...
                  ) || ( ...
                    iscell(val{1}) && (isempty(val{1}{1}) || ( ...
                      isstruct(val{1}{1}) && isempty(fieldnames(val{1}{1})) ...
                    )) ...
                  ) ...
                ) ...
              )
              error('toml:RedefinedArray', ...
                    'Arrays cannot be redefined.')
            elseif isstruct(val)
              error('toml:NameCollision', ...
                    'Table definitions cannot override existing arrays.')
            end
          otherwise
            if isstruct(val)
              error('toml:RedefinedTable', ...
                    'Tables cannot be redefined.')
            end
            error('toml:RedefinedKey', ...
                  'Keys cannot be redefined.')
        end
      end
      obj.(indx{1}) = val;
    elseif iscell(obj)
      if ischar(indx{1})
        obj{end}.(indx{1}) = val;
      else
        obj{indx{1}} = val;
      end
    end
  else
    try
      orig = toml.get_nested_field(obj, indx(1));
    catch
      if ischar(indx{2})
        orig = struct();
      else
        orig = {};
      end
    end
    new = toml.set_nested_field(orig, indx(2:end), val);
    obj = toml.set_nested_field(obj, indx(1), new);
  end
end

% GET_NESTED_FIELD get a value from inside a struct
%
%   GET_NESTED_FIELD(obj, indx) follows the pointer `indx` into the
%   struct or cell `obj`, and retrieves the value at the pointed
%   location, if it exists. If any member of the pointer chain does not
%   exist, it raises the error 'toml:NoSuchIndex'.
%
%   See also SET_NESTED_FIELD

function value = get_nested_field(obj, indx)
  % check for existence
  switch class(obj)
    case 'cell'
      if numel(obj) < indx{1}
        error('toml:NoSuchIndex', 'This index does not exist.')
      end
    case 'struct'
      if ~isfield(obj, indx{1})
        error('toml:NoSuchIndex', 'This index does not exist.')
      end
  end

  % retrieve it
  value = toml.get_item(obj, indx{1});
  % recurse if necessary
  if numel(indx) > 1
    value = toml.get_nested_field(value, indx(2:end));
  end
end

function val = get_item(obj, indx2)
  switch class(obj)
    case 'cell'
      val = obj{indx2};
    case 'struct'
      val = obj.(indx2);
  end
end


% ADJUST_KEY_STACK fix inconsistencies in a pointer into a struct
%
%   ADJUST_KEY_STACK(obj, key_stack) follows the elements of `key_stack`
%   through the structure of `obj`. If an element is missing from the
%   pointer, it is inserted in the appropriate order. The (potentially)
%   corrected pointer is returned.

function key_stack = adjust_key_stack(obj, key_stack)
  if ~isempty(key_stack)
    % see if another one is there
    switch class(obj)
      case 'cell'
        % insert missing numeric index
        if ischar(key_stack{1})
          key_stack = [{length(obj)}, key_stack];
        end
        nested_obj = obj{key_stack{1}};
      case 'struct'
        try
          nested_obj = obj.(key_stack{1});
        catch
          if numel(key_stack) > 1
            if ischar(key_stack{2})
              nested_obj = struct();
            else
              nested_obj = {};
            end
          else
            nested_obj = [];
          end
        end
    end

    % go down another level
    key_stack = [key_stack(1), toml.adjust_key_stack(nested_obj, key_stack(2:end))];
  end
end


function out = bracketarray(in)
  % BRACKETARRAY TOML representation of a MATLAB numerical array, in the style of a numpy array
  %
  %   BRACKETARRAY(array) returns the TOML/numpy-like representation of `array' with nested brackets
  %   [1, 2; 3, 4] becomes [[1,2],[3,4]]

  %   From: https://stackoverflow.com/questions/57438523/in-matlab-how-can-i-write-out-a-multidimensional-array-as-a-string-that-looks-li/57445408#57445408
  %   By: matlabbit

  out = permute(in, [2, 1, 3:ndims(in)]);
  out = string(out);

  dimsToCat = ndims(out);
  if iscolumn(out)
    dimsToCat = dimsToCat - 1;
  end

  for iDim = 1:dimsToCat
    out = "[" + join(out, ",", iDim) + "]" ;
  end
  out = char(out);
end

%% Test Suite
% disp({1, isequal(bracketarray(ones(1,1)), '[1]')})
% disp({2, isequal(bracketarray(ones(2,1)), '[[1],[1]]')})
% disp({3, isequal(bracketarray(ones(1,2)), '[1,1]')})
% disp({4, isequal(bracketarray(ones(2,2)), '[[1,1],[1,1]]')})
% disp({5, isequal(bracketarray(ones(3,2)), '[[1,1],[1,1],[1,1]]')})
% disp({6, isequal(bracketarray(ones(2,3)), '[[1,1,1],[1,1,1]]')})
% disp({7, isequal(bracketarray(ones(1,1,2)), '[[[1]],[[1]]]')})
% disp({8, isequal(bracketarray(ones(2,1,2)), '[[[1],[1]],[[1],[1]]]')})
% disp({9, isequal(bracketarray(ones(1,2,2)), '[[[1,1]],[[1,1]]]')})
% disp({10,isequal(bracketarray(ones(2,2,2)), '[[[1,1],[1,1]],[[1,1],[1,1]]]')})
% disp({11,isequal(bracketarray(ones(1,1,1,2)), '[[[[1]]],[[[1]]]]')})
% disp({12,isequal(bracketarray(ones(2,1,1,2)), '[[[[1],[1]]],[[[1],[1]]]]')})
% disp({13,isequal(bracketarray(ones(1,2,1,2)), '[[[[1,1]]],[[[1,1]]]]')})
% disp({14,isequal(bracketarray(ones(1,1,2,2)), '[[[[1]],[[1]]],[[[1]],[[1]]]]')})
% disp({15,isequal(bracketarray(ones(2,1,2,2)), '[[[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]]]')})
% disp({16,isequal(bracketarray(ones(1,2,2,2)), '[[[[1,1]],[[1,1]]],[[[1,1]],[[1,1]]]]')})
% disp({17,isequal(bracketarray(ones(2,2,2,2)), '[[[[1,1],[1,1]],[[1,1],[1,1]]],[[[1,1],[1,1]],[[1,1],[1,1]]]]')})
% disp({18,isequal(bracketarray(permute(reshape([1:16],2,2,2,2),[2,1,3,4])), '[[[[1,2],[3,4]],[[5,6],[7,8]]],[[[9,10],[11,12]],[[13,14],[15,16]]]]')})
% disp({19,isequal(bracketarray(ones(1,1,1,1,2)), '[[[[[1]]]],[[[[1]]]]]')})
% 
% assert(isequal(toml.encode(struct('a', ones(1,1))), 'a = 1'), 'ones(1,1) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(2,1))), 'a = [[1],[1]]'), 'ones(2,1) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(1,2))), 'a = [1,1]'), 'ones(1,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(2,2))), 'a = [[1,1],[1,1]]'), 'ones(2,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(3,2))), 'a = [[1,1],[1,1],[1,1]]'), 'ones(3,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(2,3))), 'a = [[1,1,1],[1,1,1]]'), 'ones(2,3) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(1,1,2))), 'a = [[[1]],[[1]]]'), 'ones(1,1,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(2,1,2))), 'a = [[[1],[1]],[[1],[1]]]'), 'ones(2,1,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(1,2,2))), 'a = [[[1,1]],[[1,1]]]'), 'ones(1,2,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(2,2,2))), 'a = [[[1,1],[1,1]],[[1,1],[1,1]]]'), 'ones(2,2,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(1,1,1,2))), 'a = [[[[1]]],[[[1]]]]'), 'ones(1,1,1,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(2,1,1,2))), 'a = [[[[1],[1]]],[[[1],[1]]]]'), 'ones(2,1,1,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(1,2,1,2))), 'a = [[[[1,1]]],[[[1,1]]]]'), 'ones(1,2,1,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(1,1,2,2))), 'a = [[[[1]],[[1]]],[[[1]],[[1]]]]'), 'ones(1,1,2,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(2,1,2,2))), 'a = [[[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]]]'), 'ones(2,1,2,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(1,2,2,2))), 'a = [[[[1,1]],[[1,1]]],[[[1,1]],[[1,1]]]]'), 'ones(1,2,2,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(2,2,2,2))), 'a = [[[[1,1],[1,1]],[[1,1],[1,1]]],[[[1,1],[1,1]],[[1,1],[1,1]]]]'), 'ones(2,2,2,2) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', permute(reshape([1:16],2,2,2,2),[2,1,3,4]))), 'a = [[[[1,2],[3,4]],[[5,6],[7,8]]],[[[9,10],[11,12]],[[13,14],[15,16]]]]'), 'permute(reshape([1:16],2,2,2,2),[2,1,3,4]) not encoded correctly!');
% assert(isequal(toml.encode(struct('a', ones(1,1,1,1,2))), 'a = [[[[[1]]]],[[[[1]]]]]'), 'ones(1,1,1,1,2) not encoded correctly!');
% 
% assert(isequal(toml.decode('a = 1'), struct('a', ones(1,1))), 'ones(1,1) not decoded correctly!');
% assert(isequal(toml.decode('a = [[1],[1]]'), struct('a', ones(2,1))), 'ones(2,1) not decoded correctly!');
% assert(isequal(toml.decode('a = [1,1]'), struct('a', ones(1,2))), 'ones(1,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[1,1],[1,1]]'), struct('a', ones(2,2))), 'ones(2,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[1,1],[1,1],[1,1]]'), struct('a', ones(3,2))), 'ones(3,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[1,1,1],[1,1,1]]'), struct('a', ones(2,3))), 'ones(2,3) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[1]],[[1]]]'), struct('a', ones(1,1,2))), 'ones(1,1,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[1],[1]],[[1],[1]]]'), struct('a', ones(2,1,2))), 'ones(2,1,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[1,1]],[[1,1]]]'), struct('a', ones(1,2,2))), 'ones(1,2,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[1,1],[1,1]],[[1,1],[1,1]]]'), struct('a', ones(2,2,2))), 'ones(2,2,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[[1]]],[[[1]]]]'), struct('a', ones(1,1,1,2))), 'ones(1,1,1,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[[1],[1]]],[[[1],[1]]]]'), struct('a', ones(2,1,1,2))), 'ones(2,1,1,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[[1,1]]],[[[1,1]]]]'), struct('a', ones(1,2,1,2))), 'ones(1,2,1,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[[1]],[[1]]],[[[1]],[[1]]]]'), struct('a', ones(1,1,2,2))), 'ones(1,1,2,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[[1],[1]],[[1],[1]]],[[[1],[1]],[[1],[1]]]]'), struct('a', ones(2,1,2,2))), 'ones(2,1,2,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[[1,1]],[[1,1]]],[[[1,1]],[[1,1]]]]'), struct('a', ones(1,2,2,2))), 'ones(1,2,2,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[[1,1],[1,1]],[[1,1],[1,1]]],[[[1,1],[1,1]],[[1,1],[1,1]]]]'), struct('a', ones(2,2,2,2))), 'ones(2,2,2,2) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[[1,2],[3,4]],[[5,6],[7,8]]],[[[9,10],[11,12]],[[13,14],[15,16]]]]'), struct('a', permute(reshape([1:16],2,2,2,2),[2,1,3,4]))), 'permute(reshape([1:16],2,2,2,2),[2,1,3,4]) not decoded correctly!');
% assert(isequal(toml.decode('a = [[[[[1]]]],[[[[1]]]]]'), struct('a', ones(1,1,1,1,2))), 'ones(1,1,1,1,2) not decoded correctly!');

% REPR TOML representation of a MATLAB object
%
%   REPR(obj) returns the TOML representation of `obj`.

function str = repr(obj, parent)

  if ispc
    newline = sprintf('\r\n');
  else
    newline = sprintf('\n');
  end

  switch class(obj)

    % strings
    case 'char'
      if isrow(obj) || isempty(obj)
        % escape characters
        esc_char = {'\\', '\b', '\t', '\n', '\f', '\r', '\"'};
        for ii = 1:length(esc_char)
          obj = strrep(obj, sprintf(esc_char{ii}), esc_char{ii});
        end
        % convert chars above \x052f to escaped unicode \u
        unicode_chars = find(double(obj) > hex2dec('52F'));
        for jj = unicode_chars
          obj = strrep(obj, obj(jj), sprintf('\\u%04s', dec2hex(double(obj(jj)))));
        end
        % wrap in basic string double quotes
        str = ['"', obj, '"'];
      else
        str = toml.repr(reshape(cellstr(obj), 1, []));
      end

    % Booleans
    case 'logical'
      reprs = {'false', 'true'};
      str = reprs{obj + 1};

    % numbers
    case 'double'
      if numel(obj) == 1
        str = lower(num2str(obj));
      else
        str = toml.bracketarray(obj);
      end

    % cell arrays
    case 'cell'
      if all(cellfun(@isstruct, obj))
        fmtter = @(a) sprintf('[[%s]]%s%s', parent, newline, repr(a));
        cel_str = cellfun(fmtter, obj, 'uniformoutput', false);
        str = strjoin(cel_str, newline);
      else
        cel_mod = cellfun(@repr, obj, 'uniformoutput', false);
        str = ['[', strjoin(cel_mod, ', '), ']'];
      end

    % structures
    case 'struct'
      obj = toml.sortfields(obj);
      fn = fieldnames(obj);
      vals = struct2cell(obj);
      str = '';
      for indx = 1:numel(vals)
        new_parent = fn{indx};
        if isstruct(vals{indx})
          if nargin > 1
            fmt_str = ['%1$s[', parent, '.%2$s]%4$s%3$s'];
            new_parent = [parent, '.', fn{indx}];
          else
            fmt_str = '%1$s[%2$s]%4$s%3$s';
          end
        elseif iscell(vals{indx}) && all(cellfun(@isstruct, vals{indx}))
          fmt_str = '%1$s%3$s%4$s';
        else
          fmt_str = '%1$s%2$s = %3$s%4$s';
        end
        str = sprintf(fmt_str, str, fn{indx}, repr(vals{indx}, new_parent), newline);
      end

    % datetime objects
    case 'datetime'
      obj.Format = 'yyyy-MM-dd''T''HH:mm:ss.SSSSSS';
      if ~isempty(obj.TimeZone)
        obj.Format = [obj.Format, 'XXX'];
      end
      str = char(obj);

    % unrecognized type
    otherwise
      error('toml:NonEncodableType', ...
            'Cannot encode type as TOML: %s', class(obj))
  end
end


% SORTFIELDS sorts fields in a struct so nesting is correct
%
%   SORTFIELDS(struct_in) returns a field-sorted struct

function struct_out = sortfields(struct_in)
  % convert to cell
  tmp = struct2cell(struct_in);
  % find which cells are structs
  sub_structs = find(cellfun(@isstruct, tmp));
  % build function to eval if a cell is a cell of all structs
  cell_of_struct = @(cell_in) iscell(cell_in) && all(cellfun(@isstruct, cell_in));
  % get index of above evaluation
  sub_cellstructs = find(cellfun(cell_of_struct, tmp));
  % list index of struct locations
  sub_structs = [sub_structs; sub_cellstructs];
  % list index order
  new_order = [setdiff(1:numel(tmp), sub_structs), sub_structs.'];
  % set new order
  struct_out = orderfields(struct_in, new_order);
end



end   %Method

end   %Classdef