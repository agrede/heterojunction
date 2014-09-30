function [x,y] = plotBands(jn, xmult)
% plotBands plots the band structure for junction(s)
%
% Copyright (C) 2014 Alex J. Grede
% GPL v3, See LICENSE.txt for details
% This function is part of heterojunction (https://github.com/agrede/heterojunction)

  x0 = -max(jn{1}.x(jn{1}.ind0,:));
  x = [];
  y = [];
  hold on;
  for k = 1:(length(jn)-1)
    [x0,tmpx,tmpy] = plotBand(x0, -xmult, jn{k}, jn{k}.ind0);
    x = [x;tmpx];
    y = [y;tmpy];
    [x0,tmpx,tmpy] = plotBand(x0, xmult, jn{k+1}, jn{k}.ind0);
    x = [x;tmpx];
    y = [y;tmpy];
  endfor
  hold off;
endfunction
