function xLift = so_lift_x_near(x, referenceX)
%SO_LIFT_X_NEAR Choose the integer x-shift nearest referenceX.
xLift = x + round(referenceX - x);
end

