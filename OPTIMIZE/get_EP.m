function EP_30day=get_EP(t_seconds,Tdata,colnumber)

% input time t_seconds in seconds
% table T


day_number=((t_seconds-mod(t_seconds,24*3600))/(24*3600))+1;

y=Tdata(:,colnumber);
t_days=1:length(y);

[mx,ix]=min(abs(t_days-day_number));

EP_mm_day=y(ix);
EP=EP_mm_day*1e-3/(24*3600);

EP_30day = movmean(EP, [29 0]);   % window = previous 29 days + today

