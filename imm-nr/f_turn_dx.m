% Jacobian of the state transition function in reentry demo.

%
% Copyright (C) 2005-2006 Simo S�rkk�
%               2007      Jouni Hartikainen
%
% This software is distributed under the GNU General Public
% Licence (version 2 or later); please refer to the file
% Licence.txt, included with the software, for details.

function df = f_turn_dx(x,param)

  dt  = param{1};
  w = x(5);
%{
  if w == 0
      coswt = 1;
      coswto = 0;
      coswtopw = 0;

      sinwt = 0;
      sinwtpw = dt;

      dsinwtpw = 0;
      dcoswtopw = -0.5*dt^2;
  else
      coswt = cos(w*dt);
      coswto = cos(w*dt)-1;
      coswtopw = coswto/w;

      sinwt = sin(w*dt);
      sinwtpw = sinwt/w;

      dsinwtpw = (w*dt*coswt - sinwt) / (w^2);
      dcoswtopw = (-w*dt*sinwt-coswto) / (w^2);
  end
  df = zeros(5,5);

  df(1,1) = 1;
  df(1,3) = sinwtpw;
  df(1,4) = coswtopw;
  df(1,5) = dsinwtpw * x(3) + dcoswtopw * x(4);

  df(2,2) = 1;
  df(2,3) = -coswtopw;
  df(2,4) = sinwtpw;
  df(2,5) = -dcoswtopw * x(3) + dsinwtpw * x(4);

  df(3,3) = coswt;
  df(3,4) = -sinwt;
  df(3,5) = -dt * sinwt * x(3) - dt * coswt * x(4);

  df(4,3) = sinwt;
  df(4,4) = coswt;
  df(4,5) = dt * coswt * x(3) - dt * sinwt * x(4);

  df(5,5) = 1;
%}
  df = zeros(5,5);

  if w == 0
    df(1,1) = 1;
    df(1,3) = dt;
    df(1,5) = -0.5*(dt^2)*x(4);

    df(2,2) = 1;
    df(2,4) = dt;
    df(2,5) = 0.5*(dt^2)*x(3);

    df(3,3) = 1;
    df(3,5) = -dt*x(4);

    df(4,4) = 1;
    df(4,5) = dt*x(3);

    df(5,5) = 1;
  else
    df(1,1) = 1;
    df(1,3) = sin(w*dt) / w;
    df(1,4) = (cos(w*dt) - 1) / w;
    df(1,5) = (cos(w*dt)*dt*x(3) / w) - ((sin(w*dt)*x(3)) / (w^2)) - ...
      (sin(w*dt)*dt*x(4) / w) - ((cos(w*dt)-1)*x(4) / (w^2));

    df(2,2) = 1;
    df(2,3) = (1-cos(w*dt)) / w;
    df(2,4) = sin(w*dt) / w;
    df(2,5) = (sin(w*dt)*dt*x(3) / w) - ((1-cos(w*dt))*x(3) / (w^2)) + ...
      (cos(w*dt)*dt*x(4) / w) - (sin(w*dt)*x(4) / (w^2));

    df(3,3) = cos(w*dt);
    df(3,4) = -sin(w*dt);
    df(3,5) = -(sin(w*dt))*dt*x(3) - (cos(w*dt))*dt*x(4);

    df(4,3) = sin(w*dt);
    df(4,4) = cos(w*dt);
    df(4,5) = cos(w*dt)*dt*x(3) - sin(w*dt)*dt*x(4);

    df(5,5) = 1;
  end

end %EOF
