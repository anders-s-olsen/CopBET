% Function to initialize results table and other small stuff
% Copyright (C) 2023 Anders Stevnhoved Olsen & Drummond E-Wen McCulloch
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see http://www.gnu.org/licenses/.
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