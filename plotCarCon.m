function xmax = plotCarCon(x0, xmult, jn, ind0)
  % plotCarCon
  %
  % Copyright (C) 2014 Alex J. Grede
  % GPL v3, See LICENSE.txt for details
  % This function is part of heterojunction (https://github.com/agrede/heterojunction)
  xmax = sign(xmult).*max(jn.x(ind0,:))-x0;
  if (sign(x0)==-1)
    x0 = xmax;
  endif
  o1 = ones(sum(jn.ind(ind0,:)),1);
  semilogy((jn.x(ind0,logical(jn.ind(ind0,:)))-x0).*xmult, ...
           abs([jn.n(logical(jn.ind(ind0,:)),:), ...
                jn.p(logical(jn.ind(ind0,:)),:), ...
                jn.NI(logical(jn.ind(ind0,:)),:)]).*1e-6);
endfunction
