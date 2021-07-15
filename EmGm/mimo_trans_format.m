function y = mimo_trans_format(x)
y = [real(x(:,1)) imag(x(:,1)) real(x(:,2)) imag(x(:,2)) real(x(:,3)) imag(x(:,3)) real(x(:,4)) imag(x(:,4))];
end