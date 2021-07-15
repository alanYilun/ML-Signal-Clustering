function y = simo_trans_format(x)
y = [real(x(:,1)) imag(x(:,1)) real(x(:,2)) imag(x(:,2))];
end