function [ x, w ] = rule15 (  )

%*****************************************************************************80
%
%% RULE15 returns the rule of degree 15.
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
    0.7749527857778351,0.9885448991378063, ...
    -.7749527857778349,-.9885448991378063, ...
    -.9070374303651182,0.9571446613308432, ...
    0.9070374303651184,-.9571446613308430, ...
    -.4303978306869286,0.9769578054468787, ...
    0.4303978306869287,-.9769578054468787, ...
    -.9756646723906326,0.1107064048513496, ...
    0.9756646723906326,-.1107064048513495, ...
    -.7388921437312957,0.7868610204187212, ...
    0.7388921437312957,-.7868610204187212, ...
    0.1995220876718269,0.6659287668239546, ...
    -.1995220876718268,-.6659287668239546, ...
    -.1934983412061240,0.8412271039808018, ...
    0.1934983412061241,-.8412271039808018, ...
    0.4882189227791580,0.8922368778153702, ...
    -.4882189227791579,-.8922368778153702, ...
    -.5772265461040059,0.9526539504944950, ...
    0.5772265461040061,-.9526539504944950, ...
    -.4474426063114782,0.5675455860909890, ...
    0.4474426063114783,-.5675455860909890, ...
    -.7044956995149931E-01,0.3256679896817100, ...
    0.7044956995149934E-01,-.3256679896817100 ];
  ys = [ ...
    -.9885448991378063,0.7749527857778350, ...
    0.9885448991378063,-.7749527857778348, ...
    -.9571446613308433,-.9070374303651183, ...
    0.9571446613308431,0.9070374303651185, ...
    -.9769578054468787,-.4303978306869286, ...
    0.9769578054468787,0.4303978306869287, ...
    -.1107064048513496,-.9756646723906326, ...
    0.1107064048513495,0.9756646723906326, ...
    -.7868610204187212,-.7388921437312957, ...
    0.7868610204187212,0.7388921437312957, ...
    -.6659287668239546,0.1995220876718268, ...
    0.6659287668239546,-.1995220876718268, ...
    -.8412271039808018,-.1934983412061240, ...
    0.8412271039808018,0.1934983412061241, ...
    -.8922368778153702,0.4882189227791580, ...
    0.8922368778153702,-.4882189227791578, ...
    -.9526539504944950,-.5772265461040060, ...
    0.9526539504944950,0.5772265461040063, ...
    -.5675455860909890,-.4474426063114783, ...
    0.5675455860909890,0.4474426063114784, ...
    -.3256679896817100,-.7044956995149933E-01, ...
    0.3256679896817100,0.7044956995149936E-01 ];
  ws = [ ...
    0.1443015463807196E-01,0.1443015463807196E-01, ...
    0.1443015463807196E-01,0.1443015463807196E-01, ...
    0.1816242330920956E-01,0.1816242330920956E-01, ...
    0.1816242330920956E-01,0.1816242330920956E-01, ...
    0.1290815898308381E-01,0.1290815898308381E-01, ...
    0.1290815898308381E-01,0.1290815898308381E-01, ...
    0.3010764365372140E-01,0.3010764365372140E-01, ...
    0.3010764365372140E-01,0.3010764365372140E-01, ...
    0.6540469907131932E-01,0.6540469907131932E-01, ...
    0.6540469907131932E-01,0.6540469907131932E-01, ...
    0.1197895531736646,0.1197895531736646, ...
    0.1197895531736646,0.1197895531736646, ...
    0.8473841548096289E-01,0.8473841548096289E-01, ...
    0.8473841548096289E-01,0.8473841548096289E-01, ...
    0.6453833756714425E-01,0.6453833756714425E-01, ...
    0.6453833756714425E-01,0.6453833756714425E-01, ...
    0.2403055376316494E-01,0.2403055376316494E-01, ...
    0.2403055376316494E-01,0.2403055376316494E-01, ...
    0.1196130510491228,0.1196130510491228, ...
    0.1196130510491228,0.1196130510491228, ...
    0.1533837904970821,0.1533837904970821, ...
    0.1533837904970821,0.1533837904970821 ];

  x(1,:) = xs;
  x(2,:) = ys;
  w = ws;

  return
end