function print_summary(header,footer,varargin)
ncols = numel(varargin);
w=12;
padding=sprintf('%s',repmat('-',1,w));
if header
    fprintf('|%s|\n|',repmat(padding,1,ncols+4));
end
for i=1:ncols
    padlen = w-numel(varargin{i});
    rempad =floor(padlen/2);

    fprintf('%s',repmat(' ',1,rempad+1));
    if mod(i,2)==0
        fprintf('%s',varargin{i});
    else
        fprintf('%s:',varargin{i});
    end
    fprintf('%s',repmat(' ',1,rempad-1));

end
fprintf('\n');
if footer
    fprintf('|%s|\n',repmat(padding,1,ncols));
end
end