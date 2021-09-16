function M = loadmatrix(hObject, eventdata, handles)

n = 6;
m = 4;
poslist = get(handles.m1, 'Position');
nlist = length(handles.m1);
M = zeros(m,n);
for k = 1:nlist,
  ks = get(handles.m1(k), 'UserData');
  M(ks(1), ks(2)) = str2num(get(handles.m1(k)), 'String');
end
