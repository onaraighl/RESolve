function [q_vec]=postprocessRainfall(t,myTable)

t_seconds=t*24*3600;
q_vec=0*t;

for i=1:length(t_seconds)
    t_val=t_seconds(i);
    q_val=get_rainfall(t_val,myTable);
    q_vec(i)=q_val;
end