% test some matlab time conversions

% restart
close all; clear; clc;

% get current time as a UNIX timestamp
% make sure to get datetime in UTC (vs local time zone)
% can verify against: https://www.unixtimestamp.com/index.php
currentDateTime = datetime('now','TimeZone','UTC');
disp(['Current datetime:      ' datestr(currentDateTime)]);
currentUnixTime = posixtime(currentDateTime);
fprintf('Current UNIX time:     %20.9f\n',currentUnixTime);  % display to nanoseconds, not sure about actual underlying accuracy

% convert UNIX time back to MATLAB datetime
% assuming UNIX time was reported in UTC
convertedDateTime = datetime(currentUnixTime,'ConvertFrom','posixtime','TimeZone','UTC');
disp(['Converted datetime:    ' datestr(convertedDateTime)]);

% get offset from UTC (i.e. GMT)
localTZOffset = tzoffset(datetime('today','TimeZone','local'));
disp(['Local offset from UTC: ' char(localTZOffset)]);