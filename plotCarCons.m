function [x,y] = plotCarCons(jn, xmult)
  % plotCarCons
  %
  % Copyright (C) 2014 Alex J. Grede
  % GPL v3, See LICENSE.txt for details
  % This function is part of heterojunction (https://github.com/agrede/heterojunction)

  x0 = -max(jn{1}.x(jn{1}.ind0,:));
  x = [];
  y = [];
  hold on;
  for k = 1:(length(jn)-1)
    [x0,tmpx,tmpy] = plotCarCon(x0, -xmult, jn{k}, jn{k}.ind0);
    x = [x;flipdim(tmpx,1)];
    y = [y;flipdim(tmpy,1)];
    [x0,tmpx,tmpy] = plotCarCon(x0, xmult, jn{k+1}, jn{k}.ind0);
    x = [x;flipdim(tmpx,1)];
    y = [y;flipdim(tmpy,1)];
  endfor
  legend({"Gamma", "Chi", "Lambda", "HH", "LH", "SO", "Impurity"})
  hold off;
endfunction
