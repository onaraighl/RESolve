function Ep=get_EP(t_seconds,metData)

Tyear=metData.Year;
EpMax=metData.EpMax;

Ep=EpMax*(cos(2*pi*t_seconds/Tyear)+1)/2; % daily EP

% EP_30day = movmean(EP, [29 0]);   % window = previous 29 days + today

