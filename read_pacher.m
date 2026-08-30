function D = read_pacher(csvfile, xlsfile)
%READ_PACHER  Parse Pacher et al. Data S1 (events) and, optionally, Data S2
%             (shoal structure). Shared by make_fig4 and make_figS4/S5/S7.
%
%  Data S1 is semicolon-delimited with 19 lines of legend and the header on
%  line 20. Seventeen columns. t_first is in FRAMES at 25 fps and missing
%  values are the literal string NA. Line 2 declares the sheet ATTACKS only,
%  which is misleading, since the 81 flybys sit in the same block with
%  risk = 'flyby'. Below the real data is a second stacked block (26 blanks,
%  a repeated header, and 47 rows shifted one column left) which must be
%  dropped; filtering on risk removes it.
%
%  Columns used: 1 file, 10 risk, 11 area (sqm), 12 true_positive,
%                13 t_first (frames), 14 bout_id.
%
%  COLUMN INDEX. The recording identifier is field 1, not field 8. Field 8 is
%  "location", of which there are only three, so keying clusters on it gave
%  31 groups (3 locations x 11 bout values) instead of 73 and understated the
%  flyby stratum at 6 rather than 26.
%
%  Returns D with fields file, attack, area, resp, tfirst (s), clust, and,
%  when xlsfile is given, D.s2 with area, N, nnd_mm, bl_mm, site.
%
%  D.tdet is the DETECTION latency vector: attacks only, NA dropped. Use it
%  for anything modeling E[T_(1) | H=1]. D.tfirst carries timings for every
%  scored event, INCLUDING 23 flybys that were answered; those are
%  false-alarm response latencies and must not enter the pooling-count fit.
%  Mixing them gives 148 latencies instead of 125.
%
%  Reproduces 177 attacks / 127 responses, 81 flybys / 26 responses,
%  125 timings, mean t_first 4.92 s.

if nargin < 1 || isempty(csvfile), csvfile = 'sciadv_adt8600_data_s1.csv'; end
HDRLINE = 20;  FPS = 25;

fid = fopen(csvfile,'r');
if fid < 0, error('read_pacher:csv','cannot open %s', csvfile); end
raw = textscan(fid,'%s','Delimiter','\n','Whitespace','');
fclose(fid);
rows = raw{1}(HDRLINE+1:end);

nm = {}; risk = {}; area = []; tp = []; tf = []; bout = {};
for i = 1:numel(rows)
    f = strsplit(rows{i}, ';', 'CollapseDelimiters', false);
    if numel(f) < 14, continue; end
    r = strtrim(f{10});
    if ~(strcmp(r,'predator_attack') || strcmp(r,'flyby')), continue; end
    nm{end+1,1}   = strtrim(f{1});               %#ok<AGROW>
    risk{end+1,1} = r;                           %#ok<AGROW>
    area(end+1,1) = na2num(f{11});               %#ok<AGROW>
    tp(end+1,1)   = na2num(f{12});               %#ok<AGROW>
    tf(end+1,1)   = na2num(f{13}) / FPS;         %#ok<AGROW>
    bout{end+1,1} = strtrim(f{14});              %#ok<AGROW>
end

D.file   = nm;
D.attack = strcmp(risk,'predator_attack');
D.area   = area;
D.resp   = tp;
D.tfirst = tf;

% bout_id RESTARTS within each recording, so it must be combined with the
% file column. Its values are '1'-'9' on attacks and 'kk' or 'other' on
% flybys. Over 18 recordings the key file|bout gives 73 clusters, 47 carrying attacks
% and 26 carrying flybys.
n = numel(nm); D.clust = cell(n,1); k = 0;
for i = 1:n
    if isempty(bout{i}) || strcmpi(bout{i},'NA')
        k = k + 1;  D.clust{i} = sprintf('%s|solo%04d', nm{i}, k);
    else
        D.clust{i} = sprintf('%s|%s', nm{i}, bout{i});
    end
end

% detection latencies: attacks only
D.tdet = D.tfirst(D.attack & ~isnan(D.tfirst));

D.s2 = [];
if nargin > 1 && ~isempty(xlsfile) && exist(xlsfile,'file')
    T = readtable(xlsfile);
    D.s2 = struct('area', T.area_sqm_, 'N', double(T.Ntags_total), ...
                  'nnd_mm', T.meanNND_mm, 'bl_mm', T.MeanBL_mm_SL_, ...
                  'site', double(T.Site_id));
end
end

function v = na2num(s)
s = strtrim(s);
if isempty(s) || strcmpi(s,'NA'), v = NaN; else, v = str2double(s); end
end
