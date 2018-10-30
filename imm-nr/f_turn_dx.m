function dF = f_turn_dx(x, param)

  dt = param{1};
  w = x(5);
  dF = zeros(5,5);

  if w < 1e-10
    dF(1,1) = 1;
    dF(1,3) = dt;
    dF(1,5) = -0.5*(dt^2)*x(4);

    dF(2,2) = 1;
    dF(2,4) = dt;
    dF(2,5) = 0.5*(dt^2)*x(3);

    dF(3,3) = 1;
    dF(3,5) = -dt*x(4);

    dF(4,4) = 1;
    dF(4,5) = dt*x(3);

    dF(5,5) = 1;
  else
    dF(1,1) = 1;
    dF(1,3) = sin(w*dt) / w;
    dF(1,4) = (cos(w*dt) - 1) / w;
    dF(1,5) = cos(w*dt)*dt*x(3) / w - sin(w*dt)*x(3) / (w^2) - ...
      sin(w*dt)*dt*x(4) / w - (cos(w*dt) - 1)*x(4) / (w^2);

    dF(2,2) = 1;
    dF(2,3) = (1 - cos(w*dt)) / w;
    dF(2,4) = sin(w*dt) / w;
    dF(2,5) = sin(w*dt)*dt*x(3) / w - (1-cos(w*dt))*x(3) / (w^2) + ...
      cos(w*dt)*dt*x(4) / w - sin(w*dt)*x(4) / (w^2);

    dF(3,3) = cos(w*dt);
    dF(3,4) = -sin(w*dt);
    dF(3,5) = -sin(w*dt)*dt*x(3) - cos(w*dt)*dt*x(4);

    dF(4,3) = sin(w*dt);
    dF(4,4) = cos(w*dt);
    dF(4,5) = cos(w*dt)*dt*x(3) - sin(w*dt)*dt*x(4);

    F(5,5) = 1;
  end

end %EOF
