function x_k = f_turn(x, param)

  dt = param{1};
  w = x(5);
  F = zeros(5,5);

  if w == 0
    F(1,1) = 1;
    F(1,3) = dt;

    F(2,1) = 1;
    F(2,4) = dt;

    F(3,3) = 1;
    F(4,4) = 1;
  else
    F(1,1) = 1;
    F(1,3) = sin(w*dt) / w;
    F(1,4) = (cos(w*dt) - 1) / w;

    F(2,2) = 1;
    F(2,3) = (1 - cos(w*dt)) / w;
    F(2,4) = sin(w*dt) / w;

    F(3,3) = cos(w*dt);
    F(3,4) = -sin(w*dt);

    F(4,3) = sin(w*dt);
    F(4,4) = cos(w*dt);

    F(5,5) = 1;
  end

  x_k = F*x;

end %EOF