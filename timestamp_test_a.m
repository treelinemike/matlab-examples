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

% show as local time
localDateTime = datetime(currentUnixTime,'ConvertFrom','posixtime','TimeZone','local');
disp(['Local datetime:        ' datestr(localDateTime)]);

% query actual local offset from UTC (i.e. GMT)
localTZOffset = tzoffset(datetime('today','TimeZone','local'));
disp(['Queried offset:        ' char(duration(localTZOffset,'Format','hh:mm:ss'))]);

% compute local offset from UTC
% to do this correctly we need to work with UNZONED datetimes (i.e. set the 'TimeZone' property to '')
localTZOffset_computed_wrong = localDateTime - convertedDateTime;  % NOTE: This will be zero because the this is the same instant in time measured at two different locations!
disp(['Wrong offset:          ' char(duration(localTZOffset_computed_wrong,'Format','hh:mm:ss'))]);
localTZOffset_computed_correct = datetime(localDateTime,'TimeZone','') - datetime(convertedDateTime,'TimeZone','');  % NOTE: This will be zero because the this is the same instant in time measured at two different locations!
disp(['Correct offset:        ' char(duration(localTZOffset_computed_correct,'Format','hh:mm:ss'))]);

% add 60 seconds to UNIX time and convert
% verifies that we can do math in seconds on unix timestamps
updatedUnixTime = currentUnixTime + 60;
updatedDateTime = datetime(updatedUnixTime,'ConvertFrom','posixtime','TimeZone','UTC');
disp(['Updated datetime:      ' datestr(updatedDateTime)]);