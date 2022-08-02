function [t,x]=ornstein_uhlenbeck_euler_maruyama ( theta, mu, sigma, x0, tmax, P, n, r)

%*****************************************************************************80
%
%% ornstein_uhlenbeck_euler_maruyama() applies Euler-Maruyama to the Ornstein-Uhlenbeck SDE.
%
%  Discussion:
%
%    The stochastic differential equation (SDE) is:
%
%      dx = theta * ( mu - x(t) ) dt + sigma dW,   
%      x(0) = x0,
%
%    The discretized Brownian path uses a constant stepsize.
%
%    A "large" time step DT_LARGE is used for the smooth portion of the
%    equation, while a smaller time step DT_SMALL is used for the
%    discretized Brownian path.  We take R small steps to equal one 
%    large step, so that:
%
%      dt_large = r * dt_small = tmax / n
%
%    For an SDE of the form:
%
%      dx = f(x(t)) dt + g(x(t)) dW(t)
%
%    the Euler-Maruyama method has the form:
%
%      x(j) = x(j-1) + f(X(j-1)) * dt_large + g(X(j-1)) * dW(j-1)
%
%    where dW(j-1) is approximated by the sum of R normal random values
%    multiplied by the square root of DT_SMALL.
%
%    Note that if SIGMA is zero, the problem becomes deterministic,
%    with solution
%
%      x(t) = mu + ( x0 - mu ) * exp ( - theta * t )
%
%  Licensing:
%
%    This code is distributed under the GNU LGPL license.
%
%  Modified:
%
%    18 January 2013
%
%  Author:
%
%    John Burkardt
%
%  Reference:
%
%    Desmond Higham,
%    An Algorithmic Introduction to Numerical Simulation of
%    Stochastic Differential Equations,
%    SIAM Review,
%    Volume 43, Number 3, September 2001, pages 525-546
%
%  Input:
%
%    real THETA, MU, SIGMA, the value of problem parameters.
%
%    real X0, the initial condition.  When studying many
%    realizations of this problem, it is usual for X0 to be chosen
%    from a normal distribution.
%
%    real TMAX, the final time.
%
%    integer N, the number of large scale time steps.
%
%    integer R, the number of small scale time steps per single
%    large scale time step.
%
%   fprintf ( 1, '\n' );
%   fprintf ( 1, 'ornstein_uhlenbeck_euler_maruyama:\n' );
%   fprintf ( 1, '  MATLAB/Octave version %s\n', version ( ) );
%   fprintf ( 1, '  Use an Euler-Maruyama method to approximate the solution of\n' );
%   fprintf ( 1, '  the Ornstein-Uhlenbeck stochastic differential equation:\n' );
%   fprintf ( 1, '\n' );
%   fprintf ( 1, '    d x(t) = theta * ( mu - x(t) ) dt + sigma dW\n' );
%   fprintf ( 1, '\n' );
%   fprintf ( 1, '  with initial condition x(0) = x0.\n' );

  if ( nargin < 1 )
    theta = 2.0;
    %fprintf ( 1, '\n' );
    %fprintf ( 1, '  Using default decay rate THETA = %g\n', theta );
  end

  if ( nargin < 2 )
    mu = 1.0;
    %fprintf ( 1, '\n' );
    %fprintf ( 1, '  Using default mean MU = %g\n', mu );
  end

  if ( nargin < 3 )
    sigma = 0.15;
    %fprintf ( 1, '\n' );
    %fprintf ( 1, '  Using default variance SIGMA = %g\n', sigma );
  end

  if ( nargin < 4 )
    x0 = 2.0;
    %fprintf ( 1, '\n' );
    %fprintf ( 1, '  Using default initial value X0 = %g\n', x0 );
  end

  if ( nargin < 5 )
    tmax = 3.0;
    %fprintf ( 1, '\n' );
    %fprintf ( 1, '  Using default final time TMAX = %g\n', tmax );
  end 

  if ( nargin < 6 )
    P= 100;
    %fprintf ( 1, '\n' );
    %fprintf ( 1, '  Using default number of large timesteps N = %d\n', n );
  end
  
  if ( nargin < 7 )
    n = 10000;
    %fprintf ( 1, '\n' );
    %fprintf ( 1, '  Using default number of large timesteps N = %d\n', n );
  end

  if ( nargin < 8 )
    r = 50;
    %fprintf ( 1, '\n' );
    %fprintf ( 1, '  Using default ratio of small to large timesteps R = %d\n', r );
  end
%
%  Set time steps.
%
  dt_large = tmax / n;
  dt_small = tmax / n / r;
%
%  Carry out the Euler-Maruyama approximate integration process.
%
  t = gpuArray.linspace ( 0, tmax, n + 1 );
  x = zeros ( P, n + 1 ,'gpuArray');

  x(:,1) = normrnd(x0,sqrt(sigma),P,1);
  for j = 1 : n
    dw = sqrt ( dt_small ) * randn ( P, r ,'gpuArray');
    x(:,j+1) = x(:,j) + dt_large * theta * ( mu - x(:,j) ) + sigma * sum ( dw(:,1:r), 2 );
  end
%
%  Plot the approximate solution.
%
  
  %xlabel ( 't', 'FontSize', 16 )
  %ylabel ( 'X(t)', 'FontSize', 16, 'Rotation', 0, 'HorizontalAlignment', 'right' )
  %title ( 'Euler-Maruyama solution of Ornstein-Uhlenbeck SDE', 'FontSize', 16 )
  %grid on

  %filename = 'ornstein_uhlenbeck_euler_maruyama.png';
  %print ( '-dpng', filename )
  %fprintf ( 1, '  Graphics saved as "%s"\n', filename );

  return
end
