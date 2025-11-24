function q=get_rainfall(t_seconds,metData)

Tyear=metData.Year;
RfMax=metData.RfMax;

q=RfMax*(cos(2*pi*t_seconds/Tyear)+1)/2; % daily Rainfall

