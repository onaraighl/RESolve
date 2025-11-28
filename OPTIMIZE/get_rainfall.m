function q=get_rainfall(t_seconds,Tdata,colnumber)

% *************************************************************************
% Matlab function get_rainfall.m
% Lennon Ó Náraigh
% November 2025
%
% Inputs:
%
% * time in seconds
% * tabular data in Tdata, where one of the cells corresponds to time in
%   DAYS, e.g. Row 1= Day 1, Row 2= Day 2, etc.
%
% Outputs:
%
% * rainfall flux q, in units of m/s.
%
% *************************************************************************

% Here, I get the day number, assuming the first entry in the table is DAY
% 1.

day_number=((t_seconds-mod(t_seconds,24*3600))/(24*3600))+1;

% Now, I extract the rainfall from the table, given that the rainfall is in
% column 6:

rf=Tdata(:,colnumber);
t_days=1:length(rf); % time in days

[mx,ix]=min(abs(t_days-day_number)); % index of array at current day

q_mm_day=rf(ix); % rainfall at current day
q=(q_mm_day*1e-3)/(24*3600); % convert into m/s.

% q_30day = movmean(q, [10 0]);   % window = previous 29 days + today
% q=q_30day;
