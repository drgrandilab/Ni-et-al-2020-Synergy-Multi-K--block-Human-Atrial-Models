%% rearrange_bars: function description
function [vec_out] = rearrange_bars(vec_in)
	vec_out = [vec_in(1:5); vec_in(11:12); vec_in(6:10); vec_in(13:end)];
