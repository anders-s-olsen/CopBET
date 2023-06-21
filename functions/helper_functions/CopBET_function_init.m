function [out,numworkers,in,NRUspecific] = CopBET_function_init(in,varargin)

parallel = true;
keepdata = true;
NRUspecific = false;

namevalue = name_value_pairs(varargin{1});
for args = 1:size(namevalue,1)
    if strcmpi(namevalue(args),'keepdata')
        keepdata = namevalue{args,2};
    elseif strcmpi(namevalue(args),'parallel')
        parallel = namevalue{args,2};
    elseif strcmpi(namevalue(args),'NRUspecific')
        NRUspecific = namevalue{args,2};
    end
    
end

if ~istable(in)
    if ismatrix(in)||ischar(in)
        % convert matrix to table with one entry
        tbl = table;
        tbl.in{1} = in;
        in = tbl;
    else
        error(['Please specify the input data as either a matrix (nxp, n>1)', ...
            'or a table of matrices tbl where the FIRST column contains the data',...
            'with a matrix for each row'])
    end
end

if parallel
    numworkers = 8;
else
    numworkers = 0;
end

if keepdata
    out = in;
else
    out = table;
end

if keepdata
    if any(strcmp(in.Properties.VariableNames,'entropy'))
        warning('Overwriting entropy column in datatbl')
    end
end

function db=name_value_pairs(c)
if isempty(c)||mod(length(c),2)==1
    db={};
    return
end
jj=1;
for ii=1:2:length(c)
     if isa(c{ii},'char')
        db{jj,1}=c{ii};db{jj,2}=c{ii+1};jj=jj+1;
     else         
         db={};
         return
     end
end