function [ file, statesRef ] = ...
    concatenateFiles(file1, file2, statesRef1, statesRef2)
%CONCATENATEFILES Concatenate file element in files.
%
%   This function will concatenate file1 and file2, and save the result in
%   file. If statesRef is provided, it will also generate the corresponding
%   state information for the files and save that into the matrix
%   statesRef.
%
%   See loadGpsData.mat for more information on the variable files and
%   collectorForStates.m for statesRef.
%
% Yaguang Zhang, Purdue, 11/28/2016

statesRef = [];

if (any(structfun(@isempty, file1)))
  for fn = fieldnames(file2)'
    file.(fn{1}) = file2.(fn{1});
  end
else
  if (strcmp(file1.type, file2.type) && strcmp(file1.id, file2.id))
      % Order the files in terms of gpsTime.
      if(file1.gpsTime(1) > file2.gpsTime(1))
          temp = file1;
          file1 = file2;
          file2 = temp;

          % For generating statesRef.
          if nargin > 2
              statesRef = [statesRef2; statesRef1];
          end
      else
          if nargin > 2
              statesRef = [statesRef1; statesRef2];
          end
      end

      % Fill the fields required.
      file.type = file1.type;
      file.id = file1.id;
      file.time = [file1.time;file2.time];
      file.gpsTime = [file1.gpsTime;file2.gpsTime];
      file.lat = [file1.lat;file2.lat];
      file.lon = [file1.lon;file2.lon];
      file.altitude = [file1.altitude;file2.altitude];
      file.speed = [file1.speed;file2.speed];
      file.bearing = [file1.bearing;file2.bearing];
      file.accuracy = [file1.accuracy;file2.accuracy];
  else
      error('Input files are not for the same vehicle!')
  end
end

end
% EOF
