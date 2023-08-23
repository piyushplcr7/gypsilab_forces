function [ x, w ] = rule12 (  )

%*****************************************************************************80
%
%% RULE12 returns the rule of degree 12.
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    07 July 2014
%
%  Author:
%
%    Original FORTRAN77 version by Hong Xiao, Zydrunas Gimbutas.
%    This MATLAB version by John Burkardt.
%
%  Reference:
%
%    Hong Xiao, Zydrunas Gimbutas,
%    A numerical algorithm for the construction of efficient quadrature
%    rules in two and higher dimensions,
%    Computers and Mathematics with Applications,
%    Volume 59, 2010, pages 663-676.
%
%  Parameters:
%
%    Input, integer N, the number of nodes.
%
%    Output, real X(2,N), the coordinates of the nodes.
%
%    Output, real W(N), the weights.
%
  xs = [ ...
    0.9711185107918885,0.6489450771045480, ...
    -.9547543104262661,-.9065777988000044, ...
    0.9288045791287373,-.9425162358139516, ...
    -.9438108523148829,-.6477285885089000, ...
    0.9399983037047607,-.6866282782659429, ...
    0.7379913501268124,-.3293152288819712, ...
    0.6556399582616308,0.4257309111871534, ...
    0.9692829476897494,-.8505721097622355, ...
    0.2079264382173936,-.8025201782903676, ...
    -.5197237466355563,-.2734035281447398E-02, ...
    -.3699658428845123,-.6558970744607242, ...
    0.3202734978144128,0.8244469498554706, ...
    0.2278925123542080E-01,0.5651896934359196E-01, ...
    0.4889968338954821,-.3555976156822369, ...
    0.7575512967066254,-.1870234315112276, ...
    -.9967741042631649 ];
  ys = [ ...
    -.8672105154213969,0.9928922644702000, ...
    0.2857493181383339,0.9656011790176721, ...
    0.8921207951072256,-.9100219543607504, ...
    0.7000393501436398,-.9664634507775797, ...
    0.4769006675305678,0.4257467094739614, ...
    -.9615153562096814,0.9693604253119810, ...
    0.7113042249747283,-.7943285461026974, ...
    -.1992840398255900,-.7990773209000775E-01, ...
    0.8876704665740045,-.6649760891057823, ...
    -.3390779542043381,-.5354471390418425, ...
    0.1502820099360215,0.8334029046137713, ...
    0.3881559105148546,-.5879856922234445, ...
    -.2000991752759776E-01,-.9646721637922943, ...
    -.2761642039851812,-.8128294162538594, ...
    0.1215344546399007,0.6390274229299343, ...
    -.4400036004541968 ];
  ws = [ ...
    0.2294319212922989E-01,0.2269718640167010E-01, ...
    0.3814000586151853E-01,0.1921567026521910E-01, ...
    0.3433859117158319E-01,0.2503589002871782E-01, ...
    0.4151906822977771E-01,0.3384747145223248E-01, ...
    0.5960510578836526E-01,0.1273691684426847, ...
    0.3629732156183973E-01,0.4288352023218015E-01, ...
    0.1186836445978463,0.1223937757234154, ...
    0.4775972626669994E-01,0.1037607311404478, ...
    0.1017344934330748,0.9441812422392200E-01, ...
    0.1662942328844954,0.1752158503094866, ...
    0.1535404788337684,0.8331450401650711E-01, ...
    0.1951461758691787,0.1055576202902579, ...
    0.1749572560557213,0.5131431669876880E-01, ...
    0.1804321296492865,0.1127530309084298, ...
    0.1400307997981144,0.1721261711453589, ...
    0.2510187133639127E-01 ];

  x(1,:) = xs;
  x(2,:) = ys;
  w = ws;

  return
end