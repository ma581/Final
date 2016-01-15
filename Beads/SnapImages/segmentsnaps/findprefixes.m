function prefixes = findprefixes;
D = dir('*-p.tif');
prefixes = [];
for i = 1:length(D),
    name = D(i).name; % e.g. D(1).name= EthoNacl001-01-p.tif
    L(i) = length(name); % e.g. 20
    p = findstr('-p.tif',name); % e.g. 15 (position where -p.tif begins)
    %       disp(name(1:p-4));
    prefixes{i} = name(1:p-4); % e.g. EthoNacl001 (cuts out -01-)
end;

delme = [];
prefixes = unique(prefixes); % finds unique prefixes

%--------------------------------------------------------------------------
%Do I really need this?
%Checks if an entry is truely unique. if not this entry is beeing deleted
%-------------------------------------------------------------------------
for i = 2:length(prefixes),
    f = strmatch(upper(prefixes(i)),upper(prefixes(1:i-1)),'exact');
    if ~isempty(f), %if not empty
        delme = [delme i];
    end;
end;


if ~isempty(delme),
    disp(['found case duplicates -- removing: ']);
    prefixes(delme),
    prefixes(delme) = [];
end;


%--------------------------------------------------------------------------
%Thos is not beeing returned out of function
%-------------------------------------------------------------------------

j = 0;

for i = 1:length(prefixes),
    D = dir([char(prefixes(i)),'-*-p.tif']);
    L(i) = length(D);
    if L(i) > 1,
        j = j + 1;
        prefixes2(j) = prefixes(i);
    end;
end;

%    prefixes = prefixes2;