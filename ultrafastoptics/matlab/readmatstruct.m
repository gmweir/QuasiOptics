function s = readmatstruct(filename)
% READMATLIST  Create structure from CSV file of materials

if nargin == 0
  filename = 'material-data.csv';
end
fid = fopen(filename);
if (fid == -1)
  error('readmatlist: Unable to open material database.');
end
fgets(fid);  % skip header line
k = 1;
count = 1;
while count
  [matname, count] = fscanf(fid, '%c', 1);
  if count ~= 0
    while matname(end) ~= ','
      matname(end + 1) = fscanf(fid, '%c', 1);
    end
    
    line = fscanf(fid, '%f,');  % read one line
    n = floor(length(line)/2);
    sell = [line(1:n,1), line([1:n] + n,1)];
    nfinite = sum(sell(:,1) ~= 0);
    
    s.(matname(1:(end-1))) = sell(1:nfinite, :);  % add field entry

    fgets(fid);  % eat up rest of line
    k = k + 1;
  end
end
fclose(fid);
